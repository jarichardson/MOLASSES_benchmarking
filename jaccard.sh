#TRUEFLOWS: western_true_flow.15m.jaccard.grd, western_true_flow.75m.jaccard.grd
#

#FILES
#in
TRUEFLOWGRD='western_true_flow.75m.jaccard.grd'
#MODLFLOWXYZ='srtm/TSRTM_4N_blindshare.xyz'
#MODLFLOWXYZ='srtm/TSRTM_4N.xyz'
#MODLFLOWXYZ='srtm/TSRTM_4N_def.xyz'
#MODLFLOWXYZ='srtm/TSRTM_4N_newparents_ps.xyz'
#MODLFLOWXYZ='srtm/TSRTM_4N_newparents_bs.xyz'
#MODLFLOWXYZ='srtm/TSRTM_4N_noparents.xyz'
#MODLFLOWXYZ='srtm/TSRTM_4N_norules.xyz'
#MODLFLOWXYZ='srtm/TSRTM_4N_norules_lim.xyz'
#MODLFLOWXYZ='srtm/TSRTM_8N.xyz'
#MODLFLOWXYZ='srtm/TSRTM_8N_def.xyz'
#MODLFLOWXYZ='srtm/TSRTM_8N_newparents_ps.xyz'
#MODLFLOWXYZ='srtm/TSRTM_8N_newparents_bs.xyz'
#MODLFLOWXYZ='srtm/TSRTM_8N_noparents.xyz'
#MODLFLOWXYZ='srtm/TSRTM_8N_norules.xyz'
MODLFLOWXYZ='tandem/TTAND_8NPARLIM_instant.xyz'



#out
MODLFLOWGRD='mdlflow.grd.tmp'
TMPSUMGRD='tmpsum.grd.tmp'
TMPINTGRD='tmpint.grd.tmp'
TMPSUMXYZ='tmpsum.xyz.tmp'


#grid flow output to same grid as the 'true flow'
echo 'Producing Modeled Flow Grid based on Observed Flow Grid'
gmtmath -C2 $MODLFLOWXYZ 0 GT = | xyz2grd -R$TRUEFLOWGRD -G$MODLFLOWGRD

#sum all locations, (1 = one hit, 2 = both hit, 0 = neither hit)
echo 'Adding Model and Observed Grids (1=one hit, 2=both hit, NaN=no hit)'
grdmath $MODLFLOWGRD $TRUEFLOWGRD ADD $TRUEFLOWGRD AND = $TMPSUMGRD

#Make all no-hits = NaN
grdclip $TMPSUMGRD -Sb0.5/NaN -G$TMPINTGRD

#output sums to ascii, suppressing all NaNs
echo 'Outputting Summed Grid to ASCII'
grd2xyz -s $TMPINTGRD > $TMPSUMXYZ

#get the number of lines (hits by either model/observation) (count)
printf '\n\nUNION of model and observed flows:\n'
wc -l $TMPSUMXYZ

#get the integrated sum (integ)
printf '\nINTEGrated Sum of union grid:\n'
gmtmath $TMPSUMXYZ SUM -o2 = | gmtinfo -Eh

printf '\nJ.I. = (integ - union) / union\n'

#number of cells in the intersection = integ - count.
#number of cells in the union = count
#JACCARD INDEX = Intesection/Union

#(integ - count) / count

./calculator_jaccard.py
echo $MODLFLOWXYZ


#RMFILES
printf '\n\nRemoving Temporary Files\n\n'
rm $MODLFLOWGRD
rm $TMPSUMGRD
rm $TMPINTGRD
rm $TMPSUMXYZ
