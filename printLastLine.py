import sys

with open(sys.argv[0]) as f:
    fd = f.fileno()
    os.lseek(fd, 0, os.SEEK_END)
    while True:
        ch = os.read(fd, 1)
        if ch == b'\n':
            line = f.read()
            break
        else:
            os.lseek(fd, -2, os.SEEK_CUR)