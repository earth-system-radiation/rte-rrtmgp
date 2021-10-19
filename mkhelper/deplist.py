#!/usr/bin/env python
import os
import re
import sys
import fnmatch
import collections

try:
    import argparse
except ImportError:
    import _argparse as argparse

_re_rule = re.compile(
    r'^[ ]*([-+\w./]+(?:[ ]+[-+\w./]+)*)[ ]*'  # targets
    r':(?:[ ]*([-+\w./]+(?:[ ]+[-+\w./]+)*))?[ ]*'  # normal prerequisites
    r'(?:\|[ ]*([-+\w./]+(?:[ ]+[-+\w./]+)*))?')  # order-only prerequisites
_meta_root = 0
_term_colors = {
    'black': 90,
    'red': 91,
    'green': 92,
    'yellow': 93,
    'blue': 94,
    'magenta': 95,
    'cyan': 96,
    'white': 97
}


def parse_args():
    class ArgumentParser(argparse.ArgumentParser):
        def convert_arg_line_to_args(self, arg_line):
            try:
                # Drop everything after the first occurrence of #:
                arg_line = arg_line[:arg_line.index('#')]
            except ValueError:
                pass

            result = []
            # Do not regard consecutive whitespaces as a single separator:
            for arg in arg_line.split(' '):
                if arg:
                    result.append(arg)
                elif result:
                    # The previous argument has a significant space:
                    result[-1] += ' '
            return result

    parser = ArgumentParser(
        fromfile_prefix_chars='@',
        description='Reads a set of MAKEFILEs and prints a topologically '
                    'sorted list of TARGETs together with their dependencies.')

    parser.add_argument(
        '-d', '--debug-file',
        help='dump debug information to DEBUG_FILE')
    parser.add_argument(
        '-t', '--target', nargs='*',
        help='names of the makefile targets; if not specified, all targets and '
             'prerequisites found in the makefiles are sent to the output')
    parser.add_argument(
        '--inc-oo', action='store_true',
        help='include order-only prerequisites in the dependency graph')
    parser.add_argument(
        '-r', '--reverse', action='store_true',
        help='print the output list in the reversed order')
    parser.add_argument(
        '--check-unique-prereq', action='append', nargs=2, metavar='PATTERN',
        help='pair of shell-like wildcards; the option enables additional '
             'consistency checks of the dependency graph: each target that '
             'matches the first pattern of the pair is checked whether it has '
             'no more than one prerequisite matching the second pattern; if '
             'the check fails, a warning message is emitted to the standard '
             'error stream')
    parser.add_argument(
        '--check-unique-basename', action='append', nargs='+',
        metavar='PATTERN',
        help='list of shell-like wildcards; the option enables additional '
             'consistency checks of the dependency graph; all targets that '
             'match at least one the patterns are checked whether none of them '
             'have the same basename; if the check fails, a warning message is '
             'emitted to the standard error stream')
    parser.add_argument(
        # Unfortunately, we cannot set nargs to 'two or more', therefore we
        # set nargs to 'one or more':
        '--check-exists-prereq', action='append', nargs='+', metavar='PATTERN',
        help='list of two or more shell-like wildcards; the option enables '
             'additional consistency checks of the dependency graph: each '
             'target that matches the first pattern of the list is checked '
             'whether it has at least one prerequisite matching any of the '
             'rest of the patterns; if the check fails, a warning message is '
             'emitted to the standard error stream')
    parser.add_argument(
        '--check-cycles', action='store_true',
        help='check whether the dependency graph is acyclic, e.g. there is no '
             'circular dependencies; if a cycle is found, a warning message is '
             'emitted to the standard output')
    parser.add_argument(
        '--check-colour', choices=_term_colors.keys(),
        help='colour the message output of the checks using ANSI escape '
             'sequences; the argument is ignored if the standard error stream '
             'is not associated with a terminal device')
    parser.add_argument(
        '-f', '--makefile', nargs='*',
        help='paths to makefiles; a single dash (-) triggers reading from '
             'the standard input stream')

    args = parser.parse_args()

    if args.check_exists_prereq:
        for pattern_list in args.check_exists_prereq:
            if len(pattern_list) < 2:
                parser.error('argument --check-exists-prereq: expected 2 or '
                             'more arguments')

    if not sys.stderr.isatty():
        args.check_colour = None

    return args


def read_makefile(makefile, inc_order_only):
    result = collections.defaultdict(list)

    if makefile == '-':
        stream = sys.stdin
    elif not os.path.isfile(makefile):
        return result
    else:
        stream = open(makefile, 'r')

    it = iter(stream)

    for line in it:
        while line.endswith('\\\n'):
            line = line[:-2]
            try:
                line += next(it)
            except StopIteration:
                break

        match = _re_rule.match(line)
        if match:
            targets = set(match.group(1).split())
            prereqs = []

            if match.group(2):
                prereqs.extend(match.group(2).split())

            if match.group(3) and inc_order_only:
                prereqs.extend(match.group(3).split())

            for target in targets:
                result[target].extend(prereqs)

    stream.close()
    return result


def visit_dfs(dep_graph, vertex,
              visited=None,
              start_visit_cb_list=None,
              finish_visit_cb_list=None,
              skip_visit_cb_list=None):
    if visited is None:
        visited = set()

    if vertex in visited:
        if skip_visit_cb_list:
            for skip_visit_cb in skip_visit_cb_list:
                skip_visit_cb(vertex)
        return

    if start_visit_cb_list:
        for start_visit_cb in start_visit_cb_list:
            start_visit_cb(vertex)

    visited.add(vertex)

    if vertex in dep_graph:
        for child in dep_graph[vertex]:
            visit_dfs(dep_graph, child, visited,
                      start_visit_cb_list,
                      finish_visit_cb_list,
                      skip_visit_cb_list)

    if finish_visit_cb_list:
        for finish_visit_cb in finish_visit_cb_list:
            finish_visit_cb(vertex)


# To make sure that the output is reproducible, we need to work with lists and
# not with sets. This helper function removes duplicates from a list while
# preserving the order.
def remove_duplicates(l):
    seen = set()
    return [x for x in l if not (x in seen or seen.add(x))]


def build_graph(makefiles, inc_oo=False):
    # Read makefiles:
    result = collections.defaultdict(list)
    for mkf in makefiles:
        mkf_dict = read_makefile(mkf, inc_oo)

        for target, prereqs in mkf_dict.items():
            result[target].extend(prereqs)

    for target in result.keys():
        result[target] = remove_duplicates(result[target])

    return result


def warn(msg, colour=None):
    sys.stderr.write("%s%s: WARNING: %s%s\n"
                     % (('\033[%dm' % _term_colors[colour]) if colour else '',
                        os.path.basename(__file__),
                        msg,
                        '\033[0m' if colour else ''))


def main():
    args = parse_args()

    if args.debug_file:
        with open(args.debug_file, 'w') as debug_file:
            debug_file.writelines([
                '# Python version: ', sys.version.replace('\n', ' '), '\n',
                '#\n',
                '# Command:\n',
                '#  ', ' '.join(sys.argv), '\n',
                '#\n',
                '# Parsed arguments:\n',
                '#  ', '\n#  '.join(
                    [k + '=' + str(v) for k, v in vars(args).items()]), '\n'])

    if args.makefile is None:
        return

    dep_graph = build_graph(args.makefile, args.inc_oo)

    if not dep_graph:
        return

    # Insert _meta_root, which will be the starting-point for the dependency
    # graph traverse:
    if args.target:
        dep_graph[_meta_root] = [t for t in args.target if t in dep_graph]
    else:
        dep_graph[_meta_root] = sorted(dep_graph.keys())

    # Visitor callbacks:
    start_visit_cb_list = []
    finish_visit_cb_list = []
    skip_visit_cb_list = []

    # Callbacks that are called once the graph is traversed:
    postprocess_cb_list = []

    if args.check_unique_prereq:
        def check_unique_prereq_start_visit_cb(vertex):
            # Skip if the vertex is _meta_root or does not have descendants:
            if vertex == _meta_root or vertex not in dep_graph:
                return
            for pattern_list in args.check_unique_prereq:
                if fnmatch.fnmatch(vertex, pattern_list[0]):
                    vertex_prereqs = dep_graph[vertex]
                    for prereq_pattern in pattern_list[1:]:
                        matching_prereqs = fnmatch.filter(vertex_prereqs,
                                                          prereq_pattern)
                        if len(matching_prereqs) > 1:
                            warn("target '%s' has more than one immediate "
                                 "prerequisite matching pattern '%s':\n\t%s"
                                 % (vertex,
                                    prereq_pattern,
                                    "\n\t".join(matching_prereqs)),
                                 args.check_colour)

        start_visit_cb_list.append(check_unique_prereq_start_visit_cb)

    if args.check_unique_basename:
        basenames = [collections.defaultdict(set) for _ in
                     range(len(args.check_unique_basename))]

        def check_unique_basename_start_visit_cb(vertex):
            # Skip if the vertex is _meta_root:
            if vertex == _meta_root:
                return
            for i, pattern_list in enumerate(args.check_unique_basename):
                for pattern in pattern_list:
                    if fnmatch.fnmatch(vertex, pattern):
                        basenames[i][os.path.basename(vertex)].add(vertex)

        start_visit_cb_list.append(check_unique_basename_start_visit_cb)

        def check_unique_basename_postprocess_cb():
            for basename_group in basenames:
                for basename, paths in basename_group.items():
                    if len(paths) > 1 and basename:
                        warn("the dependency graph contains more than one "
                             "target with basename '%s':\n\t%s"
                             % (basename, '\n\t'.join(paths)),
                             args.check_colour)

        postprocess_cb_list.append(check_unique_basename_postprocess_cb)

    if args.check_exists_prereq:
        def check_exists_prereq_start_visit_cb(vertex):
            # Skip if the vertex is _meta_root:
            if vertex == _meta_root:
                return
            for pattern_list in args.check_exists_prereq:
                if fnmatch.fnmatch(vertex, pattern_list[0]):
                    vertex_prereqs = dep_graph.get(vertex, set())
                    prereq_patterns = pattern_list[1:]
                    if not any([fnmatch.filter(vertex_prereqs, prereq_pattern)
                                for prereq_pattern in prereq_patterns]):
                        warn("target '%s' does not have an immediate "
                             "prerequisite matching any of the patterns: '%s'"
                             % (vertex,
                                "', '".join(prereq_patterns)),
                             args.check_colour)

        start_visit_cb_list.append(check_exists_prereq_start_visit_cb)

    if args.check_cycles:
        path = []

        def check_cycles_start_visit_cb(vertex):
            path.append(vertex)

        def check_cycles_skip_visit_cb(vertex):
            if vertex in path:
                start_cycle_idx = path.index(vertex)

                msg_lines = (path[1:start_cycle_idx] +
                             [path[start_cycle_idx] + ' <- start of cycle'] +
                             path[start_cycle_idx + 1:] +
                             [vertex + ' <- end of cycle'])

                warn('the dependency graph has a cycle:\n\t%s'
                     % '\n\t'.join(msg_lines), args.check_colour)

        def check_cycles_finish_visit_cb(vertex):
            path.pop()

        start_visit_cb_list.append(check_cycles_start_visit_cb)
        skip_visit_cb_list.append(check_cycles_skip_visit_cb)
        finish_visit_cb_list.append(check_cycles_finish_visit_cb)

    toposort = []

    def toposort_finish_visit_cb(vertex):
        toposort.append(vertex)
    finish_visit_cb_list.append(toposort_finish_visit_cb)

    visit_dfs(dep_graph, _meta_root,
              start_visit_cb_list=start_visit_cb_list,
              finish_visit_cb_list=finish_visit_cb_list,
              skip_visit_cb_list=skip_visit_cb_list)

    for postprocess_cb in postprocess_cb_list:
        postprocess_cb()

    # The last element of toposort is _meta_root:
    print('\n'.join(toposort[-2::-1] if args.reverse else toposort[:-1]))


if __name__ == "__main__":
    main()
