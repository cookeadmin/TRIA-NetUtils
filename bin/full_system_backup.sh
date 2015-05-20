#!/bin/bash
# sudo ./full_system_backup.sh
# When there is an error
# sent 391338506 bytes received 18231 bytes 2484804.68 bytes/sec
# total size is 1059506157988 speedup is 2707.26
# rsync warning: some files vanished before they could be transferred (code 24) at main.c(1060) [sender=3.0.7]
# or
# rsync error: some files/attrs were not transferred (see previous errors) (code 23) at main.c(1060) [sender=3.0.7]
#sent 31527578 bytes received 419 bytes 313711.41 bytes/sec
# total size is 1059506167487 speedup is 33605.25

DATE=$(date "+%s")
BACKUP_DIR=/media/cookeadmin/COOKE-LAB
BACKUP_PREFIX=$BACKUP_DIR/COOKE-LAB-BACKUP
BACKUP_CURRENT=$BACKUP_PREFIX-CURRENT
BACKUP_DATE=$BACKUP_PREFIX-$DATE
DATABASE_DIR=$BACKUP_DIR/COOKE-LAB-DATABASES
DATABASE_LIST=/TRIA-NetUtils/reference_lists/cooke-db-list

mkdir -p $BACKUP_DATE
mkdir -p $DATABASE_DIR

# RUN to perform rsync backups of entire system. Comment out if using DRY-RUN.
rsync -aAXvzm /* $BACKUP_DATE --link-dest=$BACKUP_CURRENT --delete --exclude={"/home/*/.gvfs","/home/*/.mozilla","/dev/*","/proc/*","/sys/*","/tmp/*","/run/*","/mnt/*","/media/*",/lost+found}

rm -f $BACKUP_CURRENT
ln -s $BACKUP_DATE $BACKUP_CURRENT

# BACKUP all databases using pg_dump postgresql utility.
perl /TRIA-NetUtils/bin/pg_dump_databases.pl -i $DATABASE_LIST -o $DATABASE_DIR

# DRY-RUN for testing purposes only. Comment out when you want to perform rsync backups of entire system.
# rsync -aAXvzmn --link-dest=$BACKUP_CURRENT /* $BACKUP_DATE --delete --delete-excluded --exclude={"/home/*/.gvfs","/home/*/.mozilla","/dev/*","/proc/*","/sys/*","/tmp/*","/run/*","/mnt/*","/media/*",/lost+found}

# echo "rsync completed its task without any errors or other issues." | mail -s "COOKELAB BACKUP" kevin5@ualberta.ca
