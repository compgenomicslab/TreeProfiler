
import sys
import os
import tarfile
from io import StringIO, BytesIO
import unittest
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__) + '/..'))

#from collections import namedtuple
from tempfile import NamedTemporaryFile, TemporaryDirectory

from ete4 import Tree
from treeprofiler import tree_annotate
from treeprofiler.src import utils

class TestTreeAnnotation(unittest.TestCase):
    def test_emapper(self):
        # test eggnogmapper annotation
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)

        test_tree = utils.ete4_parse("(A:1,(B:1,(E:1,D:1)Internal_1:0.5)Internal_2:0.5)Root;")
        
        # load emapper annotations
        with NamedTemporaryFile(suffix='.out.emapper.annotations') as f_annotation:
            emapper_text = '## Tue Jun  6 10:36:47 2023\n## emapper-2.1.9\n## /data/shared/home/emapper/miniconda3/envs/eggnog-mapper-2.1/bin/emapper.py --cpu 20 --mp_start_method forkserver --data_dir /dev/shm/ -o out --output_dir /emapper_web_jobs/emapper_jobs/user_data/MM_8bdu7zy0 --temp_dir /emapper_web_jobs/emapper_jobs/user_data/MM_8bdu7zy0 --override -m diamond --dmnd_ignore_warnings --dmnd_algo ctg -i /emapper_web_jobs/emapper_jobs/user_data/MM_8bdu7zy0/queries.fasta --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign denovo --num_servers 2 --report_orthologs --decorate_gff yes --excel\n##\n'
            emapper_text += '#query	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs\n'
            emapper_text += 'A	1000565.METUNv1_03972	4.99e-223	614.0	COG1348@1|root,COG1348@2|Bacteria,1MVTE@1224|Proteobacteria,2VIK4@28216|Betaproteobacteria,2KUME@206389|Rhodocyclales	206389|Rhodocyclales	P	The key enzymatic reactions in nitrogen fixation are catalyzed by the nitrogenase complex which has 2 components the iron protein and the molybdenum-iron protein	nifH	-	1.18.6.1	ko:K02588	ko00625,ko00910,ko01100,ko01120,map00625,map00910,map01100,map01120	M00175	R05185,R05496	RC00002,RC01395,RC02891	ko00000,ko00001,ko00002,ko01000	-	-	-	Fer4_NifH\n'
            emapper_text += 'B	765911.Thivi_3647	2.76e-190	530.0	COG1348@1|root,COG1348@2|Bacteria,1MVTE@1224|Proteobacteria,1RR82@1236|Gammaproteobacteria,1WW4V@135613|Chromatiales	135613|Chromatiales	P	The key enzymatic reactions in nitrogen fixation are catalyzed by the nitrogenase complex, which has 2 components the iron protein and the molybdenum-iron protein	nifH	-	1.18.6.1	ko:K02588	ko00625,ko00910,ko01100,ko01120,map00625,map00910,map01100,map01120	M00175	R05185,R05496	RC00002,RC01395,RC02891	ko00000,ko00001,ko00002,ko01000	-	-	-	Fer4_NifH\n'
            emapper_text += 'E	1009370.ALO_07448	8.77e-204	564.0	COG1348@1|root,COG1348@2|Bacteria,1TPXR@1239|Firmicutes,4H3VB@909932|Negativicutes	909932|Negativicutes	P	The key enzymatic reactions in nitrogen fixation are catalyzed by the nitrogenase complex, which has 2 components the iron protein and the molybdenum-iron protein	-	-	1.18.6.1	ko:K02588	ko00625,ko00910,ko01100,ko01120,map00625,map00910,map01100,map01120	M00175	R05185,R05496	RC00002,RC01395,RC02891	ko00000,ko00001,ko00002,ko01000	-	-	-	Fer4_NifH\n'
            emapper_text += 'D	1009370.ALO_17011	6.09e-173	483.0	COG1348@1|root,COG1348@2|Bacteria,1TPXR@1239|Firmicutes,4H2BN@909932|Negativicutes	909932|Negativicutes	P	Belongs to the NifH BchL ChlL family	-	-	1.18.6.1	ko:K02588	ko00625,ko00910,ko01100,ko01120,map00625,map00910,map01100,map01120	M00175	R05185,R05496	RC00002,RC01395,RC02891	ko00000,ko00001,ko00002,ko01000	-	-	-	Fer4_NifH\n'
            emapper_text += '## 4 queries scanned\n## Total time (seconds): 216.7490794658661\n## Rate: 9.91 q/s\n'
            f_annotation.write(emapper_text.encode())        
            f_annotation.flush()

            #metadata_dict, node_props, columns, prop2type = tree_annotate.parse_csv([f_annotation.name])
            test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
            emapper_annotations=f_annotation.name,
            )
        #print(test_tree_annotated.write(props=None,  format_root_node=True))
        
        expected_tree = Tree('(A:1[&&NHX:BRITE=ko00000|ko00001|ko00002|ko01000:BiGG_Reaction=NaN:CAZy=NaN:COG_category=P:Description=The key enzymatic reactions in nitrogen fixation are catalyzed by the nitrogenase complex which has 2 components the iron protein and the molybdenum-iron protein:EC=1.18.6.1:GOs=NaN:KEGG_Module=M00175:KEGG_Pathway=ko00625|ko00910|ko01100|ko01120|map00625|map00910|map01100|map01120:KEGG_Reaction=R05185|R05496:KEGG_TC=NaN:KEGG_ko=ko_K02588:KEGG_rclass=RC00002|RC01395|RC02891:PFAMs=Fer4_NifH:Preferred_name=nifH:eggNOG_OGs=COG1348@1|root|COG1348@2|Bacteria|1MVTE@1224|Proteobacteria|2VIK4@28216|Betaproteobacteria|2KUME@206389|Rhodocyclales:evalue=4.99e-223:max_annot_lvl=206389|Rhodocyclales:score=614.0:seed_ortholog=1000565.METUNv1_03972],(B:1[&&NHX:BRITE=ko00000|ko00001|ko00002|ko01000:BiGG_Reaction=NaN:CAZy=NaN:COG_category=P:Description=The key enzymatic reactions in nitrogen fixation are catalyzed by the nitrogenase complex_ which has 2 components the iron protein and the molybdenum-iron protein:EC=1.18.6.1:GOs=NaN:KEGG_Module=M00175:KEGG_Pathway=ko00625|ko00910|ko01100|ko01120|map00625|map00910|map01100|map01120:KEGG_Reaction=R05185|R05496:KEGG_TC=NaN:KEGG_ko=ko_K02588:KEGG_rclass=RC00002|RC01395|RC02891:PFAMs=Fer4_NifH:Preferred_name=nifH:eggNOG_OGs=COG1348@1|root|COG1348@2|Bacteria|1MVTE@1224|Proteobacteria|1RR82@1236|Gammaproteobacteria|1WW4V@135613|Chromatiales:evalue=2.76e-190:max_annot_lvl=135613|Chromatiales:score=530.0:seed_ortholog=765911.Thivi_3647],(E:1[&&NHX:BRITE=ko00000|ko00001|ko00002|ko01000:BiGG_Reaction=NaN:CAZy=NaN:COG_category=P:Description=The key enzymatic reactions in nitrogen fixation are catalyzed by the nitrogenase complex_ which has 2 components the iron protein and the molybdenum-iron protein:EC=1.18.6.1:GOs=NaN:KEGG_Module=M00175:KEGG_Pathway=ko00625|ko00910|ko01100|ko01120|map00625|map00910|map01100|map01120:KEGG_Reaction=R05185|R05496:KEGG_TC=NaN:KEGG_ko=ko_K02588:KEGG_rclass=RC00002|RC01395|RC02891:PFAMs=Fer4_NifH:Preferred_name=NaN:eggNOG_OGs=COG1348@1|root|COG1348@2|Bacteria|1TPXR@1239|Firmicutes|4H3VB@909932|Negativicutes:evalue=8.77e-204:max_annot_lvl=909932|Negativicutes:score=564.0:seed_ortholog=1009370.ALO_07448],D:1[&&NHX:BRITE=ko00000|ko00001|ko00002|ko01000:BiGG_Reaction=NaN:CAZy=NaN:COG_category=P:Description=Belongs to the NifH BchL ChlL family:EC=1.18.6.1:GOs=NaN:KEGG_Module=M00175:KEGG_Pathway=ko00625|ko00910|ko01100|ko01120|map00625|map00910|map01100|map01120:KEGG_Reaction=R05185|R05496:KEGG_TC=NaN:KEGG_ko=ko_K02588:KEGG_rclass=RC00002|RC01395|RC02891:PFAMs=Fer4_NifH:Preferred_name=NaN:eggNOG_OGs=COG1348@1|root|COG1348@2|Bacteria|1TPXR@1239|Firmicutes|4H2BN@909932|Negativicutes:evalue=6.09e-173:max_annot_lvl=909932|Negativicutes:score=483.0:seed_ortholog=1009370.ALO_17011])Internal_1:0.5[&&NHX:BRITE_counter=ko00000--2||ko00001--2||ko00002--2||ko01000--2:BiGG_Reaction_counter=NaN--2:CAZy_counter=NaN--2:COG_category_counter=P--2:Description_counter=Belongs to the NifH BchL ChlL family--1||The key enzymatic reactions in nitrogen fixation are catalyzed by the nitrogenase complex_ which has 2 components the iron protein and the molybdenum-iron protein--1:EC_counter=1.18.6.1--2:GOs_counter=NaN--2:KEGG_Module_counter=M00175--2:KEGG_Pathway_counter=ko00625--2||ko00910--2||ko01100--2||ko01120--2||map00625--2||map00910--2||map01100--2||map01120--2:KEGG_Reaction_counter=R05185--2||R05496--2:KEGG_TC_counter=NaN--2:KEGG_ko_counter=ko_K02588--2:KEGG_rclass_counter=RC00002--2||RC01395--2||RC02891--2:PFAMs_counter=Fer4_NifH--2:Preferred_name_counter=NaN--2:eggNOG_OGs_counter=1TPXR@1239|Firmicutes--2||4H2BN@909932|Negativicutes--1||4H3VB@909932|Negativicutes--1||COG1348@1|root--2||COG1348@2|Bacteria--2:evalue_avg=3.045e-173:evalue_max=6.09e-173:evalue_min=8.77e-204:evalue_std=0.0:evalue_sum=6.09e-173:max_annot_lvl_counter=909932|Negativicutes--2:score_avg=523.5:score_max=564.0:score_min=483.0:score_std=3280.5:score_sum=1047.0:seed_ortholog_counter=1009370.ALO_07448--1||1009370.ALO_17011--1])Internal_2:0.5[&&NHX:BRITE_counter=ko00000--3||ko00001--3||ko00002--3||ko01000--3:BiGG_Reaction_counter=NaN--3:CAZy_counter=NaN--3:COG_category_counter=P--3:Description_counter=Belongs to the NifH BchL ChlL family--1||The key enzymatic reactions in nitrogen fixation are catalyzed by the nitrogenase complex_ which has 2 components the iron protein and the molybdenum-iron protein--2:EC_counter=1.18.6.1--3:GOs_counter=NaN--3:KEGG_Module_counter=M00175--3:KEGG_Pathway_counter=ko00625--3||ko00910--3||ko01100--3||ko01120--3||map00625--3||map00910--3||map01100--3||map01120--3:KEGG_Reaction_counter=R05185--3||R05496--3:KEGG_TC_counter=NaN--3:KEGG_ko_counter=ko_K02588--3:KEGG_rclass_counter=RC00002--3||RC01395--3||RC02891--3:PFAMs_counter=Fer4_NifH--3:Preferred_name_counter=NaN--2||nifH--1:eggNOG_OGs_counter=1MVTE@1224|Proteobacteria--1||1RR82@1236|Gammaproteobacteria--1||1TPXR@1239|Firmicutes--2||1WW4V@135613|Chromatiales--1||4H2BN@909932|Negativicutes--1||4H3VB@909932|Negativicutes--1||COG1348@1|root--3||COG1348@2|Bacteria--3:evalue_avg=2.03e-173:evalue_max=6.09e-173:evalue_min=8.77e-204:evalue_std=0.0:evalue_sum=6.09e-173:max_annot_lvl_counter=135613|Chromatiales--1||909932|Negativicutes--2:score_avg=525.6666666666666:score_max=564.0:score_min=483.0:score_std=1654.3333333333333:score_sum=1577.0:seed_ortholog_counter=1009370.ALO_07448--1||1009370.ALO_17011--1||765911.Thivi_3647--1])Root:0[&&NHX:BRITE_counter=ko00000--4||ko00001--4||ko00002--4||ko01000--4:BiGG_Reaction_counter=NaN--4:CAZy_counter=NaN--4:COG_category_counter=P--4:Description_counter=Belongs to the NifH BchL ChlL family--1||The key enzymatic reactions in nitrogen fixation are catalyzed by the nitrogenase complex which has 2 components the iron protein and the molybdenum-iron protein--1||The key enzymatic reactions in nitrogen fixation are catalyzed by the nitrogenase complex_ which has 2 components the iron protein and the molybdenum-iron protein--2:EC_counter=1.18.6.1--4:GOs_counter=NaN--4:KEGG_Module_counter=M00175--4:KEGG_Pathway_counter=ko00625--4||ko00910--4||ko01100--4||ko01120--4||map00625--4||map00910--4||map01100--4||map01120--4:KEGG_Reaction_counter=R05185--4||R05496--4:KEGG_TC_counter=NaN--4:KEGG_ko_counter=ko_K02588--4:KEGG_rclass_counter=RC00002--4||RC01395--4||RC02891--4:PFAMs_counter=Fer4_NifH--4:Preferred_name_counter=NaN--2||nifH--2:eggNOG_OGs_counter=1MVTE@1224|Proteobacteria--2||1RR82@1236|Gammaproteobacteria--1||1TPXR@1239|Firmicutes--2||1WW4V@135613|Chromatiales--1||2KUME@206389|Rhodocyclales--1||2VIK4@28216|Betaproteobacteria--1||4H2BN@909932|Negativicutes--1||4H3VB@909932|Negativicutes--1||COG1348@1|root--4||COG1348@2|Bacteria--4:evalue_avg=1.5225e-173:evalue_max=6.09e-173:evalue_min=4.99e-223:evalue_std=0.0:evalue_sum=6.09e-173:max_annot_lvl_counter=135613|Chromatiales--1||206389|Rhodocyclales--1||909932|Negativicutes--2:score_avg=547.75:score_max=614.0:score_min=483.0:score_std=3053.5833333333335:score_sum=2191.0:seed_ortholog_counter=1000565.METUNv1_03972--1||1009370.ALO_07448--1||1009370.ALO_17011--1||765911.Thivi_3647--1];', 
                        parser=parser)   
        for leaf in test_tree_annotated.leaves():
            props = list(leaf.props.keys())
            self.assertEqual(leaf.write(props=props, parser=parser, format_root_node=True), expected_tree[leaf.name].write(props=props, parser=parser, format_root_node=True))

    def test_pfam(self):
        # test alignment
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)
        test_tree = utils.ete4_parse("(1000565.METUNv1_03972:1,(1007099.SAMN05216287:1,(1121400.SAMN02746065_101305:1,1009370.ALO_07448:1)Internal_1:0.5)Internal_2:0.5)Root;")

        # pfam data
        with NamedTemporaryFile(suffix='.out.emapper.pfam') as f_pfam:
            pfam_text = "## Tue Jun  6 10:36:50 2023\n## emapper-2.1.9\n## /data/shared/home/emapper/miniconda3/envs/eggnog-mapper-2.1/bin/emapper.py --cpu 20 --mp_start_method forkserver --data_dir /dev/shm/ -o out --output_dir /emapper_web_jobs/emapper_jobs/user_data/MM_8bdu7zy0 --temp_dir /emapper_web_jobs/emapper_jobs/user_data/MM_8bdu7zy0 --override -m diamond --dmnd_ignore_warnings --dmnd_algo ctg -i /emapper_web_jobs/emapper_jobs/user_data/MM_8bdu7zy0/queries.fasta --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign denovo --num_servers 2 --report_orthologs --decorate_gff yes --excel\n##\n# query_name	hit	evalue	sum_score	query_length	hmmfrom	hmmto	seqfrom	seqto	query_coverage\n"
            pfam_text += "1000565.METUNv1_03972	Fer4_NifH	6.3e-85	284.5	308	3	265	45	305	0.8441558441558441\n"
            pfam_text += "1007099.SAMN05216287_3993	Fer4_NifH	2.9e-138	459.5	294	1	270	4	274	0.9183673469387755\n"
            pfam_text += "1121400.SAMN02746065_101305	Oxidored_nitro	1.8e-70	237.2	734	1	398	321	721	0.5449591280653951\n"
            pfam_text += "1121400.SAMN02746065_101305	Fer4_NifH	6.5e-93	310.3	734	1	267	1	266	0.36103542234332425\n"
            pfam_text += "1009370.ALO_07448	Fer4_NifH	6.5e-98	327.1	291	1	264	6	268	0.9003436426116839\n"
            f_pfam.write(pfam_text.encode())        
            f_pfam.flush()

            with NamedTemporaryFile(suffix='.fasta') as f_msa:
                fasta_text = ">1000565.METUNv1_03972\n\
    ------------------------------------------------------------\n\
    ----------MNTVTTTHVP--LSSLKT-------------R-------R--GTS-----\n\
    QA---------D----------GEGS----VQVHQDPT--LRIGT--AKVFAVYGKGGIG\n\
    KSTTSSNLSVAFSK----L--GKRVLQIGCDPKHDSTFTLTKS-----------------\n\
    --------------------LVPTVID---ILETVD------F----------------H\n\
    S------EELRP-----EDFVFPG--YNG------------VMCVEAGG-PPAGTGCGGY\n\
    VVGQTVKLLKEHHLL--D---E--------T----DVVIFDVLGDVVCGGFAAPLQ--HA\n\
    DRALVVTANDFDSIFAMNRIVAAIQAK--SKNY---KVRLGGVIANR-----SN--AT--\n\
    -----------DQIDRFNERVGLKTMAQFPDLDV-IR--RSRL-KKATLFEM------DP\n\
    TV----E-----VEAVQHEYLRLAASLWAG----A---DPL------ECAP---MK----\n\
    DRDIFDLLGFD-------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    -------------------\n\
    >1007099.SAMN05216287_3993\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    -------------------------------------------MA--LRQCAIYGKGGIG\n\
    KSTTTQNLVSALAE----A--GQKVMIVGCDPKADSTRLILHA-----------------\n\
    -------------------KAQNSIME---MAAEAG-----------------------S\n\
    V------EDLEL-----EDVLKVG--YRD------------IKCVESGG-PEPGVGCAGR\n\
    GVITAINFLEEEGAYE-E---D--------L----DFVFYDVLGDVVCGGFAMPIRENKA\n\
    QEIYIVCSGEMMAMYAANNIAKGIVKY--ANSG---SVRLAGLICNS--RNTAR--ED--\n\
    -----------ELIMELARQLGTQMIHFVPRDNV-VQ--RAEI-RRMTVVEY--------\n\
    DP----T------AKQADEYRQLANKIVNN--R-----NFV------IPTP---IT----\n\
    MDELESLLMEFGI-----LD-----E------E---------DESII-------------\n\
    -------------------------------GKA------------AHEEAAS-A-----\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    -------------------\n\
    >1121400.SAMN02746065_101305\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------MKIAIYGKGGIG\n\
    KSTISANLSAALAK----A--GKKVLQIGCDPKHDSTRLLLGG-----------------\n\
    -------------------KRIMTALD---YMKNTP-----------------------V\n\
    G-------LQRL-----DRVLHVG--YKG------------IVCAEAGG-PEPGVGCAGR\n\
    GILSTFALFERLGLD--M---N-T------F----DVVVYDVLGDVVCGGFAVPLRQGFA\n\
    DTVYVVTSEEFMSIYAANNILKGVKNF--DQGG----HRLAGLILNS--RGTHE--NR--\n\
    -----------HPVKRFAQNVKLPVKQTVPRSEL-FR--KAEM-MEKTVVEA--------\n\
    FP----D------STEARAFHDLARDVLEN----H---TFY------PARF---LN----\n\
    EDVLEQLILQDSP-----PQ-----A------Q---------TDAEN-------------\n\
    ---------------RV------------PSKAP------------GINEKSP-EKKFKV\n\
    KSS-------------------------------------DKKSVFLSKSLLTREPLHGC\n\
    AFSGALATTTQIKDTVTVAHGPRSCTNIACQAILS-----------------AGFRLFTR\n\
    KKIL-----LE-------------------------------------------------\n\
    -------NQIAPAVISSDMDESVVIYGGKDNLVKTLEQAME-Q--------NPKAVFLVT\n\
    TCPSGVIGDDPVAAIHEIRQKYPQIPVIAVTSDGNLRGD-YMQGVLNACMEGAGALMDK-\n\
    TVTPKSHCVNILAEKNIAFNAESNFNTIADILKEMNIDI-NCRFVRNTSVEQLKGFLKAP\n\
    LNLPAYTDYFGRLMADFIDERLGIPTAKQPFPVGFSESVAWVREIADFFHES-M-AGER-\n\
    --VIDTHRRHYETMIKTYGHTLKGRRLMILTYMHNVDWIVEAAFDLG----------MEV\n\
    IKVCILNFSQD----------------NLFIT-----RYP--ERFE-VETNYDPAKRDKD\n\
    LERLKP-DLLLGNY--TPKNLPYPLHVDIIPMCPDVGFYGGLAFAHRWATLIKAPVTEGW\n\
    KNDAL--------------\n\
    >1009370.ALO_07448\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    -----------------------------------------MAKK--IKQIAIYGKGGIG\n\
    KSTTTSNISAALAV----A--GYKVMQVGCDPKSDSTNTLRGG-----------------\n\
    -------------------TYIPTVLD---TLQDRS------------------------\n\
    S--------VKL-----SEIVFEG--FHG------------VYCVEAGG-PAPGVGCAGR\n\
    GIISAVQTLKNLKVY--D---D--------L--DLDIVIYDVLGDVVCGGFAVPIREGIA\n\
    EHVFTVSSADFMAIYAANNLFKGIKKY--SNSR---GALLGGVIANS--ISAPY--AK--\n\
    -----------QIVDDFASRTKTQVVGYVPRSVT-VT--QSEL-QGKTTIEA--------\n\
    FP----D------SPQAQVYKQLAAKIAAH----E---VSA------TPSP---LE----\n\
    IEELRSWAAQWAD-----NL----VA------L---------ETGEV-------------\n\
    -------------------------------RSA------------AQSI----------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    -------------------\n"
                f_msa.write(fasta_text.encode())
                f_msa.flush()

                test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
                alignment=f_msa.name, emapper_pfam=f_pfam.name
                )

                
            expected_tree = "(1000565.METUNv1_03972:1[&&NHX:dom_arq=Fer4_NifH@171@608],(1007099.SAMN05216287:1,(1121400.SAMN02746065_101305:1[&&NHX:dom_arq=Oxidored_nitro@780@1312||Fer4_NifH@169@610],1009370.ALO_07448:1[&&NHX:dom_arq=Fer4_NifH@169@607])Internal_1:0.5[&&NHX:dom_arq=Oxidored_nitro@780@1312||Fer4_NifH@169@610])Internal_2:0.5[&&NHX:dom_arq=none@none@none])Root[&&NHX:dom_arq=Fer4_NifH@171@608];"

        self.assertEqual(test_tree_annotated.write(props=["dom_arq"], parser=parser, format_root_node=True), expected_tree)

    def test_smart(self):
        # test alignment
        # load tree
        internal_parser = "name"
        parser = utils.get_internal_parser(internal_parser)
        
        test_tree = utils.ete4_parse("(1000565.METUNv1_03972:1,(1007099.SAMN05216287:1,(1121400.SAMN02746065_101305:1,1009370.ALO_07448:1)Internal_1:0.5)Internal_2:0.5)Root;")

        # smart data
        with NamedTemporaryFile(suffix='.out.emapper.pfam') as f_smart:
            smart_text = "1000565.METUNv1_03972	AAA	41	211	0.524566042059228\n1000565.METUNv1_03972	SRP54	42	262	36307.0836782689\n1000565.METUNv1_03972	ALAD	105	217	85040.1722567821\n1000565.METUNv1_03972	LIM	135	174	1265.27135640085\n1000565.METUNv1_03972	LytTR	150	238	41327.9413786371\n1000565.METUNv1_03972	VHP	188	217	949.65191865172\n1000565.METUNv1_03972	SAF	199	257	80999.5566918776\n"
            smart_text += "1007099.SAMN05216287_3993	RAS	1	134	299.463246391791\n1007099.SAMN05216287_3993	AAA	2	130	1.12006109933941\n1007099.SAMN05216287_3993	SRP54	3	171	19746.377786172\n1007099.SAMN05216287_3993	YL1_C	4	34	62142.1549437987\n1007099.SAMN05216287_3993	UDPG_MGDP_dh_C	7	106	120144.095567571\n1007099.SAMN05216287_3993	Dak1_2	8	277	100474.2546051\n1007099.SAMN05216287_3993	HTH_DEOR	55	101	1708.14400541016\n1007099.SAMN05216287_3993	DHHA2	58	173	85384.7641607122\n1007099.SAMN05216287_3993	HhH2	78	113	701.373137135878\n1007099.SAMN05216287_3993	cNMP	105	227	547.330514787662\n1007099.SAMN05216287_3993	MGS	135	210	82904.7219178997\n1007099.SAMN05216287_3993	GATase_5	138	255	159773.057819153\n1007099.SAMN05216287_3993	DUF3585	150	270	30661.2406626041\n1007099.SAMN05216287_3993	Ribosomal_S13_N	156	217	162292.32589339\n1007099.SAMN05216287_3993	Amb_V_allergen	156	190	74760.8280411016\n1007099.SAMN05216287_3993	Malic_M	158	266	77604.2997814085\n1007099.SAMN05216287_3993	ADSL_C	167	246	61745.9686502824\n1007099.SAMN05216287_3993	NGN	177	272	881.485864401489\n"
            smart_text += "1121400.SAMN02746065_101305	AAA	2	355	8.86546657676417\n1121400.SAMN02746065_101305	H4	158	226	2188.17758527356\n1121400.SAMN02746065_101305	MeTrc	189	414	894.348474942383\n1121400.SAMN02746065_101305	FIST_C	251	453	171927.833933555\n1121400.SAMN02746065_101305	SET	294	384	1449.49197333925\n1121400.SAMN02746065_101305	GHA	299	385	559.252746615624\n1121400.SAMN02746065_101305	Cadherin_pro	350	413	90154.3519668998\n1121400.SAMN02746065_101305	IMPDH	352	522	52956.4312062364\n1121400.SAMN02746065_101305	MoCF_biosynth	353	477	15920.4774286842\n1121400.SAMN02746065_101305	ALAD	359	504	31454.7850654048\n1121400.SAMN02746065_101305	PBP5_C	362	449	10258.5684918211\n1121400.SAMN02746065_101305	MAPKK1_Int	423	511	108239.883996528\n1121400.SAMN02746065_101305	BRIGHT	543	610	514.785975446317\n"
            smart_text += "1009370.ALO_07448	AAA	4	194	16.4901572434334\n1009370.ALO_07448	SRP54	5	202	70300.2139690043\n1009370.ALO_07448	FtsA	22	194	117837.675850668\n1009370.ALO_07448	DHDPS	51	252	101873.191011028\n1009370.ALO_07448	DHHA2	99	252	43393.7684825686\n1009370.ALO_07448	ETF	103	264	43090.4551064702\n1009370.ALO_07448	GATase_5	134	291	96046.6250578487\n1009370.ALO_07448	MyTH4	139	266	525.097143529794\n1009370.ALO_07448	Haem_bd	177	273	102676.956308069\n1009370.ALO_07448	DSRM	193	254	1406.06053620549\n"
            f_smart.write(smart_text.encode())        
            f_smart.flush()

            with NamedTemporaryFile(suffix='.fasta') as f_msa:
                fasta_text = ">1000565.METUNv1_03972\n\
    ------------------------------------------------------------\n\
    ----------MNTVTTTHVP--LSSLKT-------------R-------R--GTS-----\n\
    QA---------D----------GEGS----VQVHQDPT--LRIGT--AKVFAVYGKGGIG\n\
    KSTTSSNLSVAFSK----L--GKRVLQIGCDPKHDSTFTLTKS-----------------\n\
    --------------------LVPTVID---ILETVD------F----------------H\n\
    S------EELRP-----EDFVFPG--YNG------------VMCVEAGG-PPAGTGCGGY\n\
    VVGQTVKLLKEHHLL--D---E--------T----DVVIFDVLGDVVCGGFAAPLQ--HA\n\
    DRALVVTANDFDSIFAMNRIVAAIQAK--SKNY---KVRLGGVIANR-----SN--AT--\n\
    -----------DQIDRFNERVGLKTMAQFPDLDV-IR--RSRL-KKATLFEM------DP\n\
    TV----E-----VEAVQHEYLRLAASLWAG----A---DPL------ECAP---MK----\n\
    DRDIFDLLGFD-------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    -------------------\n\
    >1007099.SAMN05216287_3993\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    -------------------------------------------MA--LRQCAIYGKGGIG\n\
    KSTTTQNLVSALAE----A--GQKVMIVGCDPKADSTRLILHA-----------------\n\
    -------------------KAQNSIME---MAAEAG-----------------------S\n\
    V------EDLEL-----EDVLKVG--YRD------------IKCVESGG-PEPGVGCAGR\n\
    GVITAINFLEEEGAYE-E---D--------L----DFVFYDVLGDVVCGGFAMPIRENKA\n\
    QEIYIVCSGEMMAMYAANNIAKGIVKY--ANSG---SVRLAGLICNS--RNTAR--ED--\n\
    -----------ELIMELARQLGTQMIHFVPRDNV-VQ--RAEI-RRMTVVEY--------\n\
    DP----T------AKQADEYRQLANKIVNN--R-----NFV------IPTP---IT----\n\
    MDELESLLMEFGI-----LD-----E------E---------DESII-------------\n\
    -------------------------------GKA------------AHEEAAS-A-----\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    -------------------\n\
    >1121400.SAMN02746065_101305\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------MKIAIYGKGGIG\n\
    KSTISANLSAALAK----A--GKKVLQIGCDPKHDSTRLLLGG-----------------\n\
    -------------------KRIMTALD---YMKNTP-----------------------V\n\
    G-------LQRL-----DRVLHVG--YKG------------IVCAEAGG-PEPGVGCAGR\n\
    GILSTFALFERLGLD--M---N-T------F----DVVVYDVLGDVVCGGFAVPLRQGFA\n\
    DTVYVVTSEEFMSIYAANNILKGVKNF--DQGG----HRLAGLILNS--RGTHE--NR--\n\
    -----------HPVKRFAQNVKLPVKQTVPRSEL-FR--KAEM-MEKTVVEA--------\n\
    FP----D------STEARAFHDLARDVLEN----H---TFY------PARF---LN----\n\
    EDVLEQLILQDSP-----PQ-----A------Q---------TDAEN-------------\n\
    ---------------RV------------PSKAP------------GINEKSP-EKKFKV\n\
    KSS-------------------------------------DKKSVFLSKSLLTREPLHGC\n\
    AFSGALATTTQIKDTVTVAHGPRSCTNIACQAILS-----------------AGFRLFTR\n\
    KKIL-----LE-------------------------------------------------\n\
    -------NQIAPAVISSDMDESVVIYGGKDNLVKTLEQAME-Q--------NPKAVFLVT\n\
    TCPSGVIGDDPVAAIHEIRQKYPQIPVIAVTSDGNLRGD-YMQGVLNACMEGAGALMDK-\n\
    TVTPKSHCVNILAEKNIAFNAESNFNTIADILKEMNIDI-NCRFVRNTSVEQLKGFLKAP\n\
    LNLPAYTDYFGRLMADFIDERLGIPTAKQPFPVGFSESVAWVREIADFFHES-M-AGER-\n\
    --VIDTHRRHYETMIKTYGHTLKGRRLMILTYMHNVDWIVEAAFDLG----------MEV\n\
    IKVCILNFSQD----------------NLFIT-----RYP--ERFE-VETNYDPAKRDKD\n\
    LERLKP-DLLLGNY--TPKNLPYPLHVDIIPMCPDVGFYGGLAFAHRWATLIKAPVTEGW\n\
    KNDAL--------------\n\
    >1009370.ALO_07448\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    -----------------------------------------MAKK--IKQIAIYGKGGIG\n\
    KSTTTSNISAALAV----A--GYKVMQVGCDPKSDSTNTLRGG-----------------\n\
    -------------------TYIPTVLD---TLQDRS------------------------\n\
    S--------VKL-----SEIVFEG--FHG------------VYCVEAGG-PAPGVGCAGR\n\
    GIISAVQTLKNLKVY--D---D--------L--DLDIVIYDVLGDVVCGGFAVPIREGIA\n\
    EHVFTVSSADFMAIYAANNLFKGIKKY--SNSR---GALLGGVIANS--ISAPY--AK--\n\
    -----------QIVDDFASRTKTQVVGYVPRSVT-VT--QSEL-QGKTTIEA--------\n\
    FP----D------SPQAQVYKQLAAKIAAH----E---VSA------TPSP---LE----\n\
    IEELRSWAAQWAD-----NL----VA------L---------ETGEV-------------\n\
    -------------------------------RSA------------AQSI----------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    ------------------------------------------------------------\n\
    -------------------\n"
                f_msa.write(fasta_text.encode())
                f_msa.flush()

                test_tree_annotated, annotated_prop2type = tree_annotate.run_tree_annotate(test_tree, 
                alignment=f_msa.name, emapper_smart=f_smart.name
                )
                
                
            expected_tree = "(1000565.METUNv1_03972:1[&&NHX:dom_arq=AAA@165@452||SRP54@168@530||ALAD@283@461||LIM@355@411||LytTR@370@502||VHP@427@461||SAF@438@525],(1007099.SAMN05216287:1,(1121400.SAMN02746065_101305:1[&&NHX:dom_arq=AAA@170@814||H4@437@532||MeTrc@478@960||FIST_C@588@999||SET@716@921||GHA@721@922||Cadherin_pro@809@959||IMPDH@811@1071||MoCF_biosynth@812@1025||ALAD@835@1052||PBP5_C@838@995||MAPKK1_Int@969@1059||BRIGHT@1092@1164],1009370.ALO_07448:1[&&NHX:dom_arq=AAA@165@478||SRP54@168@499||FtsA@185@478||DHDPS@220@575||DHHA2@359@575||ETF@363@603||GATase_5@409@710||MyTH4@414@605||Haem_bd@457@612||DSRM@477@580])Internal_1:0.5[&&NHX:dom_arq=AAA@170@814||H4@437@532||MeTrc@478@960||FIST_C@588@999||SET@716@921||GHA@721@922||Cadherin_pro@809@959||IMPDH@811@1071||MoCF_biosynth@812@1025||ALAD@835@1052||PBP5_C@838@995||MAPKK1_Int@969@1059||BRIGHT@1092@1164])Internal_2:0.5[&&NHX:dom_arq=none@none@none])Root[&&NHX:dom_arq=AAA@165@452||SRP54@168@530||ALAD@283@461||LIM@355@411||LytTR@370@502||VHP@427@461||SAF@438@525];"
        self.assertEqual(test_tree_annotated.write(props=["dom_arq"], parser=parser, format_root_node=True), expected_tree)

if __name__ == '__main__':
    unittest.main()
#pytest.main(['-v'])