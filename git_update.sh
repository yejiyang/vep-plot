#!/bin/sh

# $1, first input, describe this time commit
# determine whether $1 is empty or not
# if $1 is empty, then show the usage of this update.sh 
# else if $1 equal to sth, do the git update 
if [ -z "$1" ]; then
	echo *************usage ************
	echo usage: sh update.sh description
	echo usage: describe what is meaning of this commit.
	echo '$1' = $1
else 

## check current status, list untracked or modified files
#	git status
## check the changed stuff
#	git diff

### delete a file from git
## git rm
### mv a fileA fileB in git
## git mv
## git log
## git log -p -2 : -p means show the commit differences; 
#-2 means show the recents 2 commits
## git log -p -2

# before push, check pull first
	git pull origin master
# then do push
	git add .
	git commit -m "$1"
	git push origin master
fi
