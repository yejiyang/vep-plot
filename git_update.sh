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
#	'--all', is important, otherwise, default '--ignore-removal' will be used
	git add --all .
	git commit -m "$1"
	git push -u origin master
fi

################ for first time set up
#
####Create a new repository on the command line
#touch README.md
#git init
#git add README.md
#git commit -m "first commit"
#git remote add origin git@github.com:yejiyang/vep-plot.git
#git push -u origin master
#
####Push an existing repository from the command line
#git remote add origin git@github.com:yejiyang/vep-plot.git
#git push -u origin master
#
