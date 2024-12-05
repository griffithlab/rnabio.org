
#Inspired by instructions here:
http://steveadams.io/2016/08/03/AMI-to-Docker.html

#1. Create a new EC2 instance. Make sure to use a different AMI from the one that was originally used to create the AMI you are trying to dockerize (e.g., a different version of ubuntu) otherwise it can lead to device UUID collisions. Make sure to increase root volume to a large size so that there will be space to mount the AMI root volume and then build the docker image (e.g., 200GB total size). 

#2. Start the new EC2 instance and login (take note of the region and subregion. E.g., us-east-1d

#3. Install and configure docker
sudo snap install docker
sudo groupadd docker
sudo usermod -aG docker ${USER}

#4. Reboot the instance and then log back in
sudo reboot 

#5. Log back into the instance. 

#6. Start a screen session.

#7. Find the snapshot id of the root volume of the AMI

#8. Go to snapshots and search for that id

#9. Create a new volume from the snapshot. Make sure to choose the same region and sub-region as your new instance

#10. Attach the volume to the instance. 

#11. Find the name of the ext4 root partition of the newly attached volume
sudo lsblk -f

#It might looks something like this:

NAME         FSTYPE LABEL           UUID                                 FSAVAIL FSUSE% MOUNTPOINT
loop0                                                                          0   100% /snap/amazon-ssm-agent/7628
loop1                                                                          0   100% /snap/core18/2790
loop2                                                                          0   100% /snap/lxd/24061
loop3                                                                          0   100% /snap/snapd/20290
loop4                                                                          0   100% /snap/core20/2015
nvme0n1                                                                                 
├─nvme0n1p1  ext4   cloudimg-rootfs 1cf5eaaa-7b0d-4793-8ed8-672ca895257a      6G    21% /
├─nvme0n1p14                                                                            
└─nvme0n1p15 vfat   UEFI            E9DB-1C94                              98.3M     6% /boot/efi
nvme1n1                                                                                 
├─nvme1n1p1  ext4   cloudimg-rootfs 4f575094-453e-450d-aeed-215d8cbcbf58                
├─nvme1n1p14                                                                            
└─nvme1n1p15 vfat   UEFI            60A8-A611                                           

#12. Mount the partition
sudo mount /dev/nvme1n1p1 /mnt/

#13. Create a docker image from AMI root partition
sudo tar zcpP -C /mnt/ . | docker import - griffithlab/rnabio:0.0.3

#14. Test the image
mkdir -p ~/rnabio-workspace
docker run -v ~/rnabio-workspace:/workspace:rw --user ubuntu:ubuntu -it griffithlab/rnabio:0.0.3 /bin/bash

#15. Login to dockerhub (you will need an account with permissions to push to the repo)
docker login

#16. Push the image to dockerhub
docker push griffithlab/rnabio:0.0.3

