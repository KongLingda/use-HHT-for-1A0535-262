echo "date: $1"
echo "phase: $2"
echo ${PWD}
date=$1
phase=$2
cd ${PWD}/${date}
ftgrouppha infile=LE_combine_${phase}.pha  backfile=LE_combine_${phase}_bkg.pha outfile=LE_combine_${phase}_grp.pha grouptype=min groupscale=1000
#ftgrouppha infile=LE_combine_${phase}.pha  respfile=LE_combine_${phase}.rsp  backfile=LE_combine_${phase}_bkg.pha outfile=LE_combine_${phase}_grp.pha grouptype=optsnmin groupscale=5
ftgrouppha infile=ME_combine_${phase}.pha  backfile=ME_combine_${phase}_bkg.pha outfile=ME_combine_${phase}_grp.pha grouptype=min groupscale=1000
#ftgrouppha infile=ME_combine_${phase}.pha  respfile=ME_combine_${phase}.rsp backfile=ME_combine_${phase}_bkg.pha outfile=ME_combine_${phase}_grp.pha grouptype=optsnmin groupscale=5
ftgrouppha infile=HE_combine_${phase}.pha  backfile=HE_combine_${phase}_bkg.pha outfile=HE_combine_${phase}_grp.pha grouptype=min groupscale=5000
grppha LE_combine_${phase}_grp.pha !LE_combine_${phase}_grp.pha 'systematics 0-1534 0.01' exit
grppha ME_combine_${phase}_grp.pha !ME_combine_${phase}_grp.pha 'systematics 0-1023 0.01' exit
grppha HE_combine_${phase}_grp.pha !HE_combine_${phase}_grp.pha 'systematics 0-255 0.01' exit
cd ${PWD}
