# LONG NON-CODING GENES ANNOTATION WITH RESPECT TO THE CLOSEST PROTEIN-CODING GENE

The aim of the [FEELnc_tpLevel2gnLevelClassification.R](https://github.com/tderrien/FEELnc/blob/master/scripts/FEELnc_tpLevel2gnLevelClassification.R) is to transform *transcript-level* classification from the FEELnc_classifier module to *gene-level* classification.


The following columns tags concern the FEELnc class annotation of the LncRNA gene with respect to the nearest protein coding gene (feelLncPcg):

 - feelLncPcgClassName
 - feelLncPcgClassType
 - feelLncPcgGnId
 - feelLncPcgGnName
 - feelLncPcgGnDist

## feelLncPcgClassName: 

Abbreviation of the FEELnc classification of the LNC with respect to the closest PCG

To transfer the FEELnc information from the **transcript level** to the **gene level**, an order of importance was decided.

The class names are composed of three parts: 
 - the first part (8 letters) is composed of the main class type. 
   - For the genic classes: lncgSSex, lncgSSin, lncgASex, lncgASin. 
   - For the intergenic classes: lincDivg, lincSSup, lincSSdw, lincConv.
 - the second part (4 letters) concerns only the genic classes without subtype conflicts (see below), we add one of the three subtypes: Nested (Nest), Overlapping (Ovlp) or Containing (Cont) 
 - the third part (_n.n.n or _n.1.n) indicates that there are conflicts between annotation due to several PCGs related to the LNC locus. 

Conflicts cases are of two types: the cases in which there are more than 1 annotation relative to one unique PCG (as indicated by the `n` in the middle and `1` at the end of `n.n.1`, and the case in which there more than 1 or more annotation relative to more than 1 PCG (`n.X.n` in the feelLncPcgClassType column).    
In these cases, we prioritized the annotation in the column « feelLncPcgClassName», which gives only 1 class per gene.
 -  n.n.1 case: Genics have priority over intergenics (lncg > linc).
Among the genic, exonics have priority over intronics. 
Among exonics and intronics, the subtypes nested / containing / overlapping have the same importance. They are kept if they do not produce conflicts and are removed if there are 2 or more subtypes.
 -  n.X.n case: Same order of priority as previously (n.n.1 case). 
Concerning the intergenics, there can be annotation conflicts between the several PCGs: the LNC can be classified as lincDivg with one PCG and lincConv with another one PCG (see figure 1 for an example). We prioritize the classes as following: Divg > SS > Conv. Same-strand have priority over Conv because it could suggest an error in the modelization: the LNC could be a 5'-part or 3'part of the PCG. Between the same strand up and down (lincSSup and lincSSdw), we choose the closest. The third part of the class name of these genes is either `_n.n.n`or `_n.1.n`. 

## feelLncPcgClassType: 

gives information on three fields separated by a dot `.` (X1.X2.X3) about the classification done by FEELnc of the LNC transcript relatively to the closest PCG transcript (= LNC:PCG pair):
 - X1: number of transcripts of the LNC gene: `1` if 1 transcript, `n` if more than one transcript,
 - X2: number of feelnc class(es) associated to the LNC:PCG pair: `1` if 1 class, `n` if more than one class (the `unclassified` class does not count),
 - X3: number of PCG gene(s) concerned by this (these) annotation(s): `1` if 1 PCG gene, `n` if more than one PCG gene.
1.1.1: the LNC has 1 transcript, with 1 annotation associated to 1 PCG.
n.1.1: the LNC has several transcripts, all with the same annotation associated to the same PCG.
n.n.1: the LNC has several transcripts, with different annotations associated to the same PCG.
n.n.n: the LNC has several transcripts, with different annotations associated to different PCG.
n.1.n: the LNC has several transcripts, all with the same annotation but associated to different PCG.
unclassified: the LNC is either alone in a contig (beginning by AADN. or KQ), or no interactions were found within the 100 000 pb sliding window used by FEELnc. 


![Image of FEELnc_class_conflict](http://tools.genouest.org/data/tderrien/cnrs_umr6290/FEELnc_classGn.png)

Figure. Examples of configurations corresponding to a `n.n.n` type (top) or a `n.1.n` type (bottom)

## feelLncPcgGnId: 

Unique identifier of the protein-coding gene relatively to which a LNC gene is classified by FEELnc

## feelLncPcgGnName: 

Name of the protein-coding gene relatively to which a LNC gene is classified by FEELnc

## feelLncPcgGnDist: 

Distance (in bp), as calculated by FEELnc, between the protein-coding gene relatively to which a LNC gene is classified by FEELnc and the LNC gene. 


# Contributions

Many thanks to Frederic Jehl, [Kevin Muret](https://github.com/kevmuret) and Sandrine Lagarrigue. 

# Warnings

It has been tested on chicken lncRNAs annotated by FEELnc. 
