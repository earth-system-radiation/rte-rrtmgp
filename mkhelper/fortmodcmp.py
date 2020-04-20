#!/usr/bin/env python
import sys

BUF_MAX_SIZE = 512


def skip_sequence(stream, sequence):
    """
    Finds the first occurrence of a sequence of bytes in a binary stream and
    sets the streams's current position right after it. Returns True if the
    sequence is found and False otherwise, The length of the sequence must not
    exceed BUF_MAX_SIZE.
    """
    sequence_size = len(sequence)
    while 1:
        buf = stream.read(BUF_MAX_SIZE)
        idx = buf.find(sequence)
        if idx < 0:
            if len(buf) < BUF_MAX_SIZE:
                return False
            else:
                stream.seek(1 - sequence_size, 1)
        else:
            stream.seek(idx + sequence_size - len(buf), 1)
            return True


def mods_differ(filename1, filename2, compiler_name=None):
    """
    Checks whether two Fortran module files are essentially different. Some
    compiler-specific logic is required for compilers that generate different
    module files for the same source file (e.g. the module files might contain
    timestamps). This implementation is inspired by CMake.
    """
    with open(filename1, 'rb') as stream1:
        with open(filename2, 'rb') as stream2:
            if compiler_name == "intel":
                # The first byte encodes the version of the module file format:
                if stream1.read(1) != stream2.read(1):
                    return True
                # The block before the following magic sequence contains might
                # change from compilation to compilation, probably due to a
                # second resolution timestamp in it:
                magic_sequence = b'\x0A\x00'  # the same as \n\0
                if not (skip_sequence(stream1, magic_sequence) and
                        skip_sequence(stream2, magic_sequence)):
                    return True
            elif compiler_name == "gnu":
                magic_sequence = b'\x1F\x8b'  # the magic number of gzip
                if stream1.read(2) == magic_sequence:  # gfortran 4.9 or later
                    stream1.seek(0)
                else:  # gfortran 4.8 or older
                    # Skip the first line in the text file containing a
                    # timestamp:
                    magic_sequence = b'\x0A'  # the same as \n
                    if not (skip_sequence(stream1, magic_sequence) and
                            skip_sequence(stream2, magic_sequence)):
                        return True
            elif compiler_name == "portland":
                for _ in range(2):  # the first two lines must be identical
                    if stream1.readline(BUF_MAX_SIZE) != \
                            stream2.readline(BUF_MAX_SIZE):
                        return True
                # The next line is a timestamp followed by the sequence
                # '\nenduse\n':
                magic_sequence = b'\x0A\x65\x6E\x64\x75\x73\x65\x0A'
                if not (skip_sequence(stream1, magic_sequence) and
                        skip_sequence(stream2, magic_sequence)):
                    return True
            # Compare the rest (or everything for unknown compilers):
            while 1:
                buf1 = stream1.read(BUF_MAX_SIZE)
                buf2 = stream2.read(BUF_MAX_SIZE)
                if buf1 != buf2:
                    return True
                if not buf1:
                    return False


# We try to make this as fast as possible, therefore we do not parse arguments
# properly:
exit(mods_differ(sys.argv[1], sys.argv[2],
                 sys.argv[3].lower() if len(sys.argv) > 3 else None))
