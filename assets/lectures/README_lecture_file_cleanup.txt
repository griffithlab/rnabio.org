Note, periodically it may be necessary to remove and clean up large files to keep git repo manageable.
This was done Nov 2023 using guidance from here:
https://www.phase2technology.com/blog/removing-large-files-git-bfg
https://www.tempertemper.net/blog/changing-your-git-history

#Clone or fetch/pull the repo to make sure you are up to date
cd ~/git
git clone git@github.com:griffithlab/rnabio.org.git

#First, copy the old files to genomedata
cd rnabio.org/assets/lectures/cshl
tar -zcvf 2022.tar.gz 2022/
scp -i [genomedata.pem] 2022.tar.gz ubuntu@[ip_address]:/data/rnaseq-tutorial/lectures/cshl/

#Log in to genomedata.org and unpack the lecture files and make sure in appropriate place
ssh -i [genomedata.pem] ubuntu@[ip_address]
cd /data/rnaseq-tutorial/lectures/cshl/
tar -zxvf 2022.tar.gz
exit

#Delete the files from the repo (note - this won't actually reduce the repo size since the files remain in git history)
cd ~/git/rnabio.org
git rm -r assets/lectures/cshl/2022/
git commit -m 'delete legacy files'
git push

#Take note of size of repo and then delete for now
cd ~/git/
du -h rnabio.org
rm -rf rnabio.org

#Install bfg
brew install bfg

#Clone the repo as a mirror - this will create a copy of the .git directory instead of usual set of files
git clone --mirror git@github.com:griffithlab/rnabio.org.git

#Create a backup of this repo
cp -r rnabio.org.git rnabio.org.git-backup

#Remove files from git history - note you can only specify folder names, not paths. Beware of duplicate folder names
bfg --delete-folders 2022 rnabio.org.git

#Expire and prune the git history 
git reflog expire --expire=now --all && git gc --prune=now --aggressive

#Push changes to remote - Note, you may see a lot of "remote rejected" warnings. This is not necessarily a problem
git push

#Go to github, navigate back in history to some recent commits, before the files were deleted, to confirm they are gone from history

#Get rid of mirror version of repo
cd ~/git
rm -rf rnabio.org.git

#Re-clone the repo as normal, check size for a reduction
cd ~/git
git clone git@github.com:griffithlab/rnabio.org.git
du -h rnabio.org
