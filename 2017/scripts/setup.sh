# Mount HPC home to ~/hpc via sshfs
cd; mkdir -p hpc; fusermount -u hpc
sshfs -o follow_symlinks login2.hpc.qmul.ac.uk:/data/home/$USER hpc
