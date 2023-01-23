## HemeMap Data

The methodology was described in [our paper](https://www.nature.com/articles/s41590-022-01370-4).
 
All interactions across specturm of hematopoiesis can be downloaded [HERE](https://osf.io/nkcdq).
In this Rdata, there are two objects that are row-matched, `interation_df2` is the interaction table without filtering (the last 3 columns describing coordinates of cis element)
`interaction_geoATACmean18` is the activity score for each interaction across 18 heme cell types. If you just want all the interactions specifically activated in HSC, you can simply select all the lines that interaction_geoATACmean18$HSC>8.91

If you have any questions, please reach out to Dr. Fulong Yu (fyu@broadinstitute.org).
