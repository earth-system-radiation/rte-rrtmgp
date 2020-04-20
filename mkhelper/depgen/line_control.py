import os
import re

from depgen import file_in_dir


class LCProcessor:
    _re_lc = re.compile(r'^#\s*[1-9]\d*\s*"(.*?)"\s*(?:[1-9]\d*)?')

    def __init__(self, stream,
                 include_roots=None):
        self.include_roots = include_roots

        # Callbacks:
        self.lc_callback = None
        self.debug_callback = None

        self._stream = stream

    def readline(self):
        while 1:
            line = self._stream.readline()
            if not line:
                return line

            match = LCProcessor._re_lc.match(line)
            if match:
                filepath = match.group(1)
                if os.path.isfile(filepath):
                    if not self.include_roots or any(
                            [file_in_dir(filepath, d)
                             for d in self.include_roots]):
                        if self.lc_callback:
                            self.lc_callback(filepath)
                        if self.debug_callback:
                            self.debug_callback(
                                line, 'accepted file \'%s\'' % filepath)
                    elif self.debug_callback:
                        self.debug_callback(
                            line, 'ignored (file \'%s\' '
                                  'is not in the source roots)' % filepath)
                elif self.debug_callback:
                    self.debug_callback(line, 'ignored (file not found)')
                continue

            return line

    @property
    def name(self):
        return self._stream.name

    def close(self):
        self._stream.close()
