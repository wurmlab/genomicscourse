# Fix the sshfs stuff.
set -x

cd; fusermount -u hpc; rmdir hpc