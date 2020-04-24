import re

from depgen import IncludeFinder, StreamStack, file_in_dir, open23, \
    find_unquoted_string


class Preprocessor:
    _re_ifdef = re.compile(r'^#\s*if(n)?def\s+([a-zA-Z_]\w*)')
    _re_if_expr = re.compile(r'^#\s*if((?:\s|\().*)')

    _re_elif = re.compile(r'^#\s*elif((?:\s|\().*)')
    _re_else = re.compile(r'^#\s*else(?:\s.*)')
    _re_endif = re.compile(r'^#\s*endif(?:\s.*)')

    _re_include = re.compile(r'^#\s*include\s+(?:"(.*?)"|<(.*?)>)')
    _re_define = re.compile(r'^#\s*define\s+([a-zA-Z_]\w*)(\(.*\))?\s+(.*)$')
    _re_undef = re.compile(r'^#\s*undef\s+([a-zA-Z_]\w*)')

    # matches "defined MACRO_NAME" and "defined (MACRO_NAME)"
    _re_defined_call = re.compile(
        r'(defined\s*(\(\s*)?([a-zA-Z_]\w*)(?(2)\s*\)))')

    _re_identifier = re.compile(
        r'(([a-zA-Z_]\w*)(\s*\(\s*(\w+(?:\s*,\s*\w+)*)?\s*\))?)')

    def __init__(self,
                 stream,
                 include_order=None,
                 include_sys_order=None,
                 include_dirs=None,
                 include_roots=None,
                 try_eval_expr=False,
                 inc_sys=False,
                 predefined_macros=None):
        self.include_roots = include_roots
        self.try_eval_expr = try_eval_expr
        self.inc_sys = inc_sys

        # Callbacks:
        self.include_callback = None
        self.debug_callback = None

        self._include_finder = IncludeFinder(include_order, include_dirs)
        self._include_sys_finder = IncludeFinder(include_sys_order,
                                                 include_dirs)
        self._include_stack = StreamStack()
        self._include_stack.add(stream)

        if predefined_macros:
            self._macros = dict(predefined_macros)
        else:
            self._macros = dict()

        # Stack of #if-#else blocks holds one of the following:
        # 1 - keep current branch and ignore another
        # -1 - ignore current branch and keep another
        # 0 - keep both branches (if failed to evaluate expression)
        self._if_state_stack = []

        # Each #elif is interpreted as a combination of #else and #if.
        # Thus, each #elif increments the number of #if blocks that are
        # closed with #endif statement. This numbers are stored in a separate
        # stack:
        self._states_per_endif_stack = []

    def readline(self):
        while 1:
            line = self._include_stack.readline()
            line = self._replace_continuation(line)
            line = self._remove_block_comments(line)

            if line.isspace():
                continue

            # if(n)def directive
            match = Preprocessor._re_ifdef.match(line)
            if match:
                macro, negate, state = match.group(2), bool(match.group(1)), 0
                if not self._branch_is_dead():
                    state = 1 if bool(macro in self._macros) ^ negate else -1
                    if self.debug_callback:
                        self.debug_callback(
                            line, 'evaluated to ' +
                                  ('True' if state > 0 else 'False'))
                elif self.debug_callback:
                    self.debug_callback(
                        line, 'was not evaluated (dead branch)')
                self._if(state)
                continue

            # if directive
            match = Preprocessor._re_if_expr.match(line)
            if match:
                expr, state = match.group(1), 0
                if not self._branch_is_dead():
                    if self.try_eval_expr:
                        state = self._evaluate_expr_to_state(expr)
                        if self.debug_callback:
                            self.debug_callback(
                                line, 'evaluated to ' +
                                      ('True' if state > 0 else
                                       ('False' if state < 0 else
                                        'Unknown (evaluation failed)')))
                    elif self.debug_callback:
                        self.debug_callback(
                            line, 'was not evaluated (evaluation disabled)')
                elif self.debug_callback:
                    self.debug_callback(
                        line, 'was not evaluated (dead branch)')
                self._if(state)
                continue

            # elif directive
            match = Preprocessor._re_elif.match(line)
            if match:
                self._else()
                expr, state = match.group(1), 0
                if not self._branch_is_dead():
                    if self.try_eval_expr:
                        state = self._evaluate_expr_to_state(expr)
                        if self.debug_callback:
                            self.debug_callback(
                                line, 'evaluated to ' +
                                      ('True' if state > 0 else
                                       ('False' if state < 0 else
                                        'Unknown (evaluation failed)')))
                    elif self.debug_callback:
                        self.debug_callback(
                            line, 'was not evaluated (evaluation disabled)')
                elif self.debug_callback:
                    self.debug_callback(
                        line, 'was not evaluated (dead branch)')
                self._elif(state)
                continue

            # else directive
            match = Preprocessor._re_else.match(line)
            if match:
                self._else()
                continue

            # endif directive
            match = Preprocessor._re_endif.match(line)
            if match:
                self._endif()
                continue

            if self._branch_is_dead() and self.debug_callback is None:
                continue

            # define directive
            match = Preprocessor._re_define.match(line)
            if match:
                if not self._branch_is_dead():
                    self._define(*match.group(1, 2, 3))
                    if self.debug_callback:
                        self.debug_callback(line, 'accepted')
                elif self.debug_callback:
                    self.debug_callback(line, 'ignored (dead branch)')
                continue

            # undef directive
            match = Preprocessor._re_undef.match(line)
            if match:
                if not self._branch_is_dead():
                    self._macros.pop(match.group(1), None)
                    if self.debug_callback:
                        self.debug_callback(line, 'accepted')
                elif self.debug_callback:
                    self.debug_callback(line, 'ignored (dead branch)')
                continue

            # include directive
            match = Preprocessor._re_include.match(line)
            if match:
                if not self._branch_is_dead():

                    if match.lastindex == 1:  # quoted form
                        filepath = self._include_finder.find(
                            match.group(1),
                            self._include_stack.root_name,
                            self._include_stack.current_name)
                    elif match.lastindex == 2:  # angle-bracket form
                        if self.inc_sys:
                            filepath = self._include_sys_finder.find(
                                match.group(2),
                                self._include_stack.root_name,
                                self._include_stack.current_name)
                        else:
                            if self.debug_callback:
                                self.debug_callback(line,
                                                    'ignored (system header)')
                            continue
                    else:
                        if self.debug_callback:
                            self.debug_callback(line,
                                                'ignored (internal error)')
                        continue

                    if filepath:
                        if not self.include_roots or any(
                                [file_in_dir(filepath, d)
                                 for d in self.include_roots]):
                            self._include_stack.add(open23(filepath, 'r'))
                            if self.include_callback:
                                self.include_callback(filepath)
                            if self.debug_callback:
                                self.debug_callback(
                                    line, 'included file \'%s\'' % filepath)
                        elif self.debug_callback:
                            self.debug_callback(
                                line,
                                'ignored (file \'%s\' '
                                'is not in the source roots)' % filepath)
                    elif self.debug_callback:
                        self.debug_callback(line, 'ignored (file not found)')
                elif self.debug_callback:
                    self.debug_callback(line, 'ignored (dead branch)')
                continue

            if self._branch_is_dead():
                continue

            return line

    @property
    def name(self):
        return self._include_stack.root_name

    def close(self):
        self._include_stack.clear()

    def _define(self, name, args=None, body=None):
        if name != 'defined':
            self._macros[name] = (args, '' if body is None else body)

    def _if(self, state):
        self._if_state_stack.append(state)
        self._states_per_endif_stack.append(1)

    def _else(self):
        if self._if_state_stack:
            self._if_state_stack[-1] = -self._if_state_stack[-1]

    def _elif(self, state):
        self._if_state_stack.append(state)
        self._states_per_endif_stack[-1] += 1

    def _endif(self):
        if self._if_state_stack:
            pop_count = self._states_per_endif_stack.pop()
            for _ in range(pop_count):
                self._if_state_stack.pop()

    def _branch_is_dead(self):
        return any(state < 0 for state in self._if_state_stack)

    def _replace_continuation(self, line):
        while line.endswith('\\\n'):
            suffix = self._include_stack.readline()
            line = line[:-2] + suffix
        return line

    def _remove_block_comments(self, line):
        while 1:
            # Check whether the line contains an unquoted block comment
            # initiator '/*':
            start_idx = find_unquoted_string('/*', line)
            if start_idx < 0:
                return line
            # Check whether the line contains a block comment
            # terminator '*/' (even if it is quoted):
            term_idx = line.find('*/', start_idx + 2)
            while term_idx < 0:
                # The block is not terminated yet, read the next line:
                next_line = self._include_stack.readline()
                line_length = len(line)
                if next_line:
                    line += next_line
                    term_idx = line.find('*/', line_length)
                else:
                    term_idx = line_length
            else:
                # Replace the block of comments with a single
                # space:
                line = '%s %s' % (line[:start_idx], line[term_idx + 2:])

    def _evaluate_expr_to_state(self, expr):
        prev_expr = expr
        while 1:
            # replace calls to function "defined"
            defined_calls = re.findall(Preprocessor._re_defined_call, expr)
            for call in defined_calls:
                expr = expr.replace(
                    call[0], '1' if call[2] in self._macros else '0')

            identifiers = re.findall(Preprocessor._re_identifier, expr)

            for identifier in identifiers:
                if identifier[1] == 'defined':
                    return 0

                macro = self._macros.get(identifier[1], None)
                if identifier[2]:
                    # potential call to a function
                    if macro is None:
                        # call to undefined function
                        return 0
                    elif macro[0] is not None:
                        # we can't evaluate function-like macros
                        return 0
                    else:
                        # identifier is defined as object-like macro
                        expr = expr.replace(identifier[0],
                                            macro[1] + identifier[2])
                else:
                    # no function call
                    if macro is None or macro[0] is not None:
                        # macro is not defined or
                        # defined as function-like macro
                        expr = expr.replace(identifier[0], '0')
                    else:
                        # identifier is defined as object-like macro
                        expr = expr.replace(identifier[0], macro[1])

            if prev_expr == expr:
                break
            else:
                prev_expr = expr

        expr = expr.replace('||', ' or ')
        expr = expr.replace('&&', ' and ')
        expr = expr.replace('!', 'not ')

        try:
            result = bool(eval(expr, {}))
            return 1 if result else -1
        except:
            return 0

