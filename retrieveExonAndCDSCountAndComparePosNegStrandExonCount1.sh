#! /bin/bash/
echo "Please supply  a valid file name as input"
## taking a file name as input
read fname   
if [ -z "$fname" ]
then 
     echo "Please enter a file name:"
     exit
fi    
terminal=`tty` 
exec < $fname
## initialising variables that will be used for counting the various parameters to zero
no_of_exons=0
no_of_CDS=0
not_exon_cds=0
no_of_genes=0
garbage_var=0
positive_strand_exon=0
positive_strand_not_exon=0
negative_strand_exon=0
negative_strand_not_exon=0
gene_word="uc"
## while loop to read in each column and then do to carry out each comparison
while read col1 col2 col3 col4 col5 col6 col7 col8 col9
do
    if [ "${col3}" == "exon" ]; then
    let no_of_exons="no_of_exons + 1"
    elif [ "${col3}" == "CDS" ]; then
    let no_of_CDS="$no_of_CDS + 1"
    else [ "${col3}" != "exon" ] && [ "${col3}" == "CDS" ]
     let not_exon_cds="$not_exon_cds + 1"
     fi ## end of if statement

     if [ "${col7}" == "+" ];then
       if [ "${col3}" == "exon" ];then
          let positive_strand_exon="positive_strand_exon + 1"
       fi
     else if [ "${col7}" == "-" ];then 
       if [ "${col3}" == "exon" ];then
          let negative_strand_exon="positive_strand_exon + 1"
       fi
    fi 
   fi 
done
 if [ $positive_strand_exon -gt $negative_strand_exon ]; then
    echo "The positive strand has more exons than negative strand"
 elif [ $positive_strand_exon -lt $negative_strand_exon ]; then
   echo "The negative strand has more exons than positive strand"
 else
   echo "The positive and negative strand has equal number of exons in this dataset"  
   echo "The number of exons on positive strand and negative strand are $positive_strand_exon and $negative_strand_exon respectively"  
  fi   


while read col1 col2 col3 col4 col5 col6 col7 col8 col9 
do            
        if ["${col7}" != "+" ] && ["${col7}" != "-" ];then     
            echo "Strand designation may be missing" 
        fi  
done            
                             
echo "The number of exons is $no_of_exons"
echo "The number of coding sequences or CDS is $no_of_CDS"   
echo "The number of start and stop codon  coordinates are $not_exon_cds"  
echo "The number of genes is $no_of_genes" 

 
exec < $terminal

