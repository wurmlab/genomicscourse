# Mount HPC home to ~/hpc via sshfs and make a symlink to course data directory
# in their HPC home directory (i.e., ~/hpc on the work PC).
set -x
cd; mkdir -p hpc; fusermount -u hpc
sshfs -o follow_symlinks login2.hpc.qmul.ac.uk:/data/home/$USER hpc
ssh login.hpc.qmul.ac.uk ln -fs /data/SBCS-MSc-BioInf/data/ 2017-09-BIO721_genome_bioinformatics_input
