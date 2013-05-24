#!/bin/bash

# By Ed Wiget
# This fixes dropbox sync issues on linux

mkdir -p $HOME/tmp
# get a list of files executable now
find ~/Dropbox -type f -perm -u+x > $HOME/tmp/dropbox_files-`date +%Y%m%d`

# fix the permissions
chown -R $USER ~/Dropbox
chmod -R u+rw ~/Dropbox
chown -R $USER ~/.dropbox
chmod -R u+rw ~/.dropbox
# remove any conflicting files from the file list above step 1
grep -v -e "conflicted copy" -e "Case Conflict" $HOME/tmp/dropbox_files-`date +%Y%m%d` > $HOME/tmp/dropbox_files-`date+%Y%m%d`.txt
# set the executable permissions back
for files in `echo $HOME/tmp/dropbox_files-\`date +%Y%m%d\`.txt` ; do chmod u+x "${files}" ; done

# remove any files that are in conflict
find ~/Dropbox -type f -name \*"conflicted copy"\* -exec rm -f "{}" \;
find ~/Dropbox -type f -name \*"Case Conflict"\* -exec rm -f "{}" \;

# remove temp files
rm -f $HOME/tmp/dropbox_files-`date +%Y%m%d`
rm -f $HOME/tmp/dropbox_files-`date +%Y%m%d`.txt
