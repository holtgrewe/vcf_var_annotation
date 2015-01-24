.. role:: optional
.. role:: preferred

Variant annotations in VCF format
=================================

:Version: 1.0: January 2015
:Authors: - Pablo Cingolani
          - Fiona Cunningham
          - Will McLaren
          - Kai Wang

Functional annotations field names and meanings for VCF files.
This document is intended as a standard way of representing variant annotations within VCF files (INFO fields).
It also specifies how to handle some inconsistencies, border cases and how to improve agreement with HGVS notation as well as Sequence Ontology terms and their putative impact.
The aim of this standard is to: i) make pipeline development easier, ii) facilitate benchmarking, and iii) improve some problems in edge cases

Color guide:

* **Optional** :optional:`items are highlighted in green`
* **Preferred** :preferred:`items are highlighted in yellow`
* **Mandatory** items are not highlighted

General Guidelines
------------------

We use the name "effect" and "consequence" interchangeably, meaning "functional annotation".

* VCF INFO field name ANN, stands for "annotations"

* Data fields are encoded separated by pipe sign "``|``" the order of fields is written in the VCF header.

* When comparing genomic coordinates, the comparison should be done first by chromosome names (compared alphabetically), then by start position, and finally by end position.

* Special characters: Comma, space, tab, newline or pipe characters ("``,``", "`` ``"``, "``\t``", "``\n``", "``|``", etc.) can be either:

   * .. container:: preferred

         Convert to underscore ("``_``").

     This is the preferred way.

     How about the "``p.=``" to describe synonymous variants?
     Since "``=``" is an illegal character in VCF specification, we can use an alternative notation, such as "``p.(Leu54Leu)``"

     HGVS says:

        .. container:: preferred

            Description of so called "silent" changes in the format p.(Leu54Leu) (or p.(L54L)) **should** not be used.
            When desired such changes can be described using p.(=)

     HGVS recommendation discourages the use of the format "p.(Leu54Leu)", but does not forbid it (the spec. says "should not" instead of "must not").

   * Encoded as %XX (same as URL encoding).
     This may be needed to express HGVS "p.(=)"

* Multiple "effects / consequences" are separated by comma.

   .. container:: optional

      * :optional:`Optional: Annotations are sorted by:`

           1. Effect/Consequence: Estimated deleteriousness.
              Compare using "most deleterious" when multiple consequences are predicted.

           2. In case of coding consequence: Best transcript support level (TSL http://www.ensembl.org/Help/Glossary?id=492) or Canonical transcript should be first.

           3. Feature genomic coordinates.

           4. Feature ID (compared alphabetically, even if the ID is a number).

Header fields
-------------

The variant annotations are encoded in an ``INFO`` field of name ``ANN``.
This field is defined in the header with its ``Number`` property set to ``.``, its ``Type`` set to ``String``.
The ``Description`` is a string formed as follows.::

   "Functional annotations: '${list of fields}' "

The following shows a full example of a VCF ``ANN`` header line.
The example includes line breaks that are not permitted in VCF, though.::

   ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">

The space before and after single quotation marks (``'``) can be missing , as can the space left and right of the pipe symbols (``|``).
That is, the following header line is equivalent to the earlier one.::

   ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations:'Allele |Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p|cDNA.pos / cDNA.length|CDS.pos / CDS.length|AA.pos / AA.length|Distance|ERRORS / WARNINGS / INFO'">

Field order and meaning
-----------------------

* ``Allele`` (or ``ALT``):

   * In case of multiple ALT fields, this helps to identify which ALT we are referring to.
     E.g.::

        #CHROM  POS     ID  REF ALT QUAL    FILTER  INFO
        chr1    123456  .   C   A   .       .       ANN=A|...
        chr1    234567  .   A   G,T .       .       ANN=G|..., T|...

   * In case of cancer sample, when comparing somatic versus germline using a non-standard reference (e.g. one of the ALTs is the reference) the format should be ALT-REFERENCe.
     E.g.::

        #CHROM  POS     ID  REF ALT QUAL    FILTER  INFO
        chr1    123456  .   A   C,G .       .       ANN=G-C|...

   * Compound variants: two or more variants affecting the annotations (e.g. two consecutive SNPs conforming a MNP, two consecutive frame_shift variants that "recover" the frame).
     In this case, the Allele field should include a reference to the other variant/s included in the annotation::

        #chrom  pos     id  ref alt qual    filter  info
        chr1    123456  .   A   T   .       .       ANN=|...
        chr1    123457  .   C   G   .       .       ANN=C-chr1:123456_A>T|...

* ``Annotation`` (a.k.a. effect or consequence):
  Annotated using Sequence Ontology terms.
  Multiple effects can be concatenated using "``&``".::

    #chrom  pos     id  ref alt qual    filter  info
    chr1    123456  .   a   t   .       .       ANN=A|intron_variant&nc_transcript_variant

* ``Annotation_impact``:
  A simple estimation of putative impact / deleteriousness : {HIGH, MODERATE, LOW, MODIFIER}

* ``Gene_Name``:
  Common gene name (HGNC).
  Optional: use closest gene when the variant is "intergenic".

* ``Gene_ID``:
  Gene ID

* ``Feature_type``:
  Which type of feature is in the next field (e.g. transcript, motif, miRNA, etc.).
  It :preferred:`is preferred to use Sequence Ontology (SO) terms`, but "custom" (user defined) are allowed.::

      ANN=A|stop_gained|HIGH|||transcript|...

   Tissue specific features may include cell type / tissue information separated by semicolon e.g.::

      ANN=A|histone_binding_site|LOW|||H3K4me3:HeLa-S3|...

* ``Feature_ID``:
  Depending on the annotation, this may be:
  Transcript ID (:preferred:`preferably using version number`), Motif ID, miRNA, ChipSeq peak, Histone mark, etc.
  Note: Some features may not have ID (e.g. histone marks from custom Chip范eq experiments may not have a unique ID).

* ``Transcript_BioType``:
  The bare minimum is at least a description on whether the transcript is ``{"Coding","Noncoding"}``.
  Whenever possible, use ENSEMBL biotypes.

* ``Rank``:
  Exon or Intron rank / total number of exons or introns.

* ``HGVS.c``:
  Variant using HGVS notation (DNA level)

* ``HGVS.p``:
  If variant is coding, this field describes the variant using HGVS notation (Protein level).
  Since transcript ID is already mentioned in "feature ID", it may be omitted here.

* ``cDNA_position`` :optional:`/ (cDNA_len optional)`:
  Position in cDNA and trancript's cDNA length (one based).

* ``CDS_position`` :optional:`/ (CDS_len optional)`:
  Position and number of coding bases (one based includes START and STOP codons).

* ``AA.pos`` / (Aa.length optional):
  Position and number of AA (one based, including START, but not STOP)

.. container:: optional

     * ``Distance`` to feature: All items in this field are options, so the field could be empty.

         * Up/Downstream:
           Distance to first / last codon

         * Intergenic:
           Distance to closest gene

         * Distance to closest Intron boundary in exon (+/- up/downstream).
           If same, use positive number.

         * Distance to closest exon boundary in Intron (+/- up/downstream)

         * Distance to first base in MOTIF

         * Distance to first base in miRNA

         * Distance to exon虹ntron boundary in splice_site or splice_region

         * ChipSeq peak: Distance to summit (or peak center)

         * Histone mark / Histone state: Distance to summit (or peak center)

* ``ERRORS / WARNINGS / INFO``:
  Add errors, warnings or informative message that can affect annotation accuracy.
  It can be added using either "codes" (as shown in column 1, e.g.  W1) or "message types" (as shown in column 2, e.g.  WARNING_REF_DOES_NOT_MATCH_GENOME).
  All these errors, warnings or information messages messages are optional.

  ==== ========================================= ===================
  Code Message type                              Description / Notes
  ==== ========================================= ===================
  E1   ERROR_CHROMOSOME_NOT_FOUND                Chromosome does not exists in reference genome database.
                                                 Typically indicates a mismatch between the chromosome names in the input file and the chromosome names used in the reference genome.

  E2   ERROR_OUT_OF_CHROMOSOME_RANGE             The variant's genomic coordinate is greater than chromosome's length.

  W1   WARNING_REF_DOES_NOT_MATCH_GENOME         This means that the "REF" field in the input VCF file does not match the reference genome.
                                                 This warning may indicate a conflict between input data and data from reference genome (for instance is the input VCF was aligned to a different reference genome).

  W2   WARNING_SEQUENCE_NOT_AVAILABLE            Reference sequence is not available, thus no inference could be performed.

  W3   WARNING_TRANSCRIPT_INCOMPLETE             A protein coding transcript having a non衫ultiple of 3 length.
                                                 It indicates that the reference genome has missing information about this particular transcript.

  W4   WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS   A protein coding transcript has two or more STOP codons in the middle of the coding sequence (CDS).
                                                 This should not happen and it usually means the reference genome may have an error in this transcript.

  W5   WARNING_TRANSCRIPT_NO_START_CODON         A protein coding transcript does not have a proper START codon.
                                                 It is rare that a real transcript does not have a START codon, so this probably indicates an error or missing information in the reference genome.

  I1   INFO_REALIGN_3_PRIME                      Variant has been realigned to the most 3計rime position within the transcript.
                                                 This is usually done to to comply with HGVS specification to always report the most 3計rime annotation.

  I2   INFO_COMPOUND_ANNOTATION                  This effect is a result of combining more than one variants (e.g. two consecutive SNPs that conform an MNP, or two consecutive frame_shift variants that compensate frame).

  I3   INFO_NON_REFERENCE_ANNOTATION             An alternative reference sequence was used to calculate this annotation (e.g. cancer sample comparing somatic vs. germline).
  ==== ========================================= ===================

Consistency between HGVS and functional annotations
---------------------------------------------------

In some cases there might be inconsistent reporting between 'annotation' and HGVS.
This is due to the fact that VCF recommends aligning to the leftmost coordinate, whereas HGSV recommends aligning to the "most 3計rime coordinate".

For instance, an InDel on the edge of an exon, which has an "intronic" annotation according to VCF alignment recommendation, can lead to a "stop_gained" when aligned using HGVS's recommendation (using the most 3計rime possible alignment).
So the "annotation" sub苯ield will report "intron" whereas HGVS sub苯ield will report a "stop_gained".
This is obviously inconsistent and must be avoided.

In order to report annotations that are consistent with HGVS notation, variants must be re苔ligned according to each transcript's strand (i.e. align the variant according to the transcript's most 3計rime coordinate).
Then annotations are calculated, thus the reported annotations will be consistent with HGVS notation.
:optional:`Annotation software should have a command line option to override this behaviour` (e.g. ``--no_shift_hgvs``).

A full example
--------------

To help the further understanding of the VCF field proposed here, the following shows output of SnpEff 4.1 on the example VCF provided with this tool.

::

    ##SnpEffVersion="4.1a (build 2015-01-14), by Pablo Cingolani"
    ##SnpEffCmd="SnpEff  GRCh37.75 examples/test.chr22.vcf "
    ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
    ##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' ">
    ##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' ">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    22	17071756	.	T	C	.	.	ANN=C|3_prime_UTR_variant|MODIFIER|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.*11A>G|||||11|,C|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|transcript|ENST00000430910|processed_pseudogene||n.*397A>G|||||4223|
    22	17072035	.	C	T	.	.	ANN=T|missense_variant|MODERATE|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.1406G>A|p.Gly469Glu|1666/2034|1406/1674|469/557||,T|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|transcript|ENST00000430910|processed_pseudogene||n.*397G>A|||||3944|
    22	17072258	.	C	A	.	.	ANN=A|missense_variant|MODERATE|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.1183G>T|p.Gly395Cys|1443/2034|1183/1674|395/557||,A|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|transcript|ENST00000430910|processed_pseudogene||n.*397G>T|||||3721|
    22	17072674	.	G	A	.	.	ANN=A|missense_variant|MODERATE|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.767C>T|p.Pro256Leu|1027/2034|767/1674|256/557||,A|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|transcript|ENST00000430910|processed_pseudogene||n.*397C>T|||||3305|
    22	17072747	.	T	C	.	.	ANN=C|missense_variant|MODERATE|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.694A>G|p.Met232Val|954/2034|694/1674|232/557||,C|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|transcript|ENST00000430910|processed_pseudogene||n.*397A>G|||||3232|
    22	17072781	.	C	T	.	.	ANN=T|synonymous_variant|LOW|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.660G>A|p.Pro220Pro|920/2034|660/1674|220/557||,T|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|transcript|ENST00000430910|processed_pseudogene||n.*397G>A|||||3198|
    22	17073043	.	C	T	.	.	ANN=T|missense_variant|MODERATE|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.398G>A|p.Arg133Gln|658/2034|398/1674|133/557||,T|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|transcript|ENST00000430910|processed_pseudogene||n.*397G>A|||||2936|
    22	17073066	.	A	G	.	.	ANN=G|synonymous_variant|LOW|CCT8L2|ENSG00000198445|transcript|ENST00000359963|protein_coding|1/1|c.375T>C|p.Ala125Ala|635/2034|375/1674|125/557||,G|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|transcript|ENST00000430910|processed_pseudogene||n.*397T>C|||||2913|

Annotations and putative impacts
--------------------------------

The following list describes the suggested putative impact for some Sequence Ontology terms often used in functional annotations.

=============== ======================
Putative Impact Sequence Ontology term
=============== ======================
HIGH            `chromosome_number_variation <http://sequenceontology.org/browser/current_svn/term/SO:1000182>`_
HIGH            `exon_loss_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001572>`_
HIGH            `frameshift_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001589>`_
HIGH            `rare_amino_acid_variant <http://sequenceontology.org/browser/current_svn/term/SO:0002008>`_
HIGH            `splice_acceptor_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001574>`_
HIGH            `splice_donor_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001575>`_
HIGH            `start_lost <http://sequenceontology.org/browser/current_svn/term/SO:0002012>`_
HIGH            `stop_gained <http://sequenceontology.org/browser/current_svn/term/SO:0001587>`_
HIGH            `stop_lost <http://sequenceontology.org/browser/current_svn/term/SO:0001578>`_
HIGH            `transcript_ablation <http://sequenceontology.org/browser/current_svn/term/SO:0001893>`_
MODERATE        3_prime_UTR_truncation + exon_loss
MODERATE        5_prime_UTR_truncation + exon_loss_variant
MODERATE        `coding_sequence_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001580>`_
MODERATE        `disruptive_inframe_deletion <http://sequenceontology.org/browser/current_svn/term/SO:0001826>`_
MODERATE        `disruptive_inframe_insertion <http://sequenceontology.org/browser/current_svn/term/SO:0001824>`_
MODERATE        `inframe_deletion <http://sequenceontology.org/browser/current_svn/term/SO:0001824>`_
MODERATE        `inframe_insertion <http://sequenceontology.org/browser/current_svn/term/SO:0001821>`_
MODERATE        `missense_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001583>`_
MODERATE        `regulatory_region_ablation <http://sequenceontology.org/browser/current_svn/term/SO:0001894>`_
MODERATE        `splice_region_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001630>`_
MODERATE        `TFBS_ablation <http://sequenceontology.org/browser/current_svn/term/SO:0001895>`_
LOW             5_prime_UTR_premature start_codon_gain_variant
LOW             `initiator_codon_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001582>`_
LOW             `splice_region_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001630>`_
LOW             `splice_region_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001630>`_
LOW             `start_retained <http://sequenceontology.org/browser/current_svn/term/SO:0002019>`_
LOW             `stop_retained_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001567>`_
LOW             `stop_retained_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001567>`_
LOW             `synonymous_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001819>`_
MODIFIER        `3_prime_UTR_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001624>`_
MODIFIER        `5_prime_UTR_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001623>`_
MODIFIER        `coding_sequence_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001580>`_
MODIFIER        `conserved_intergenic_variant <http://sequenceontology.org/browser/current_svn/term/SO:0002017>`_
MODIFIER        `conserved_intron_variant <http://sequenceontology.org/browser/current_svn/term/SO:0002018>`_
MODIFIER        `downstream_gene_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001632>`_
MODIFIER        `exon_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001791>`_
MODIFIER        `feature_elongation <http://sequenceontology.org/browser/current_svn/term/SO:0001907>`_
MODIFIER        `feature_truncation <http://sequenceontology.org/browser/current_svn/term/SO:0001906>`_
MODIFIER        `gene_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001564>`_
MODIFIER        `intergenic_region <http://sequenceontology.org/browser/current_svn/term/SO:0000605>`_
MODIFIER        `intragenic_variant <http://sequenceontology.org/browser/current_svn/term/SO:0002011>`_
MODIFIER        `intron_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001627>`_
MODIFIER        `mature_miRNA_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001620>`_
MODIFIER        `miRNA <http://sequenceontology.org/browser/current_svn/term/SO:0000276>`_
MODIFIER        `NMD_transcript_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001621>`_
MODIFIER        `non_coding_transcript_exon_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001792>`_
MODIFIER        `non_coding_transcript_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001619>`_
MODIFIER        `regulatory_region_amplification <http://sequenceontology.org/browser/current_svn/term/SO:0001891>`_
MODIFIER        `regulatory_region_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001566>`_
MODIFIER        `TF_binding_site_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001782>`_
MODIFIER        `TFBS_amplification <http://sequenceontology.org/browser/current_svn/term/SO:0001892>`_
MODIFIER        `transcript_amplification <http://sequenceontology.org/browser/current_svn/term/SO:0001889>`_
MODIFIER        `transcript_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001576>`_
MODIFIER        `upstream_gene_variant <http://sequenceontology.org/browser/current_svn/term/SO:0001631>`_
=============== ======================

Annotations sort order
----------------------

When comparing two annotations, the "most deleterious" one is shown first.
It is recommended annotation programs clearly state their respective "deleteriousness" order.
:preferred:`This is an example of such putative sorting order:`

#. chromosome_number_variation
#. exon_loss_variant
#. frameshift_variant
#. stop_gained
#. stop_lost
#. start_lost
#. splice_acceptor_variant
#. splice_donor_variant
#. rare_amino_acid_variant
#. missense_variant
#. inframe_insertion12. disruptive_inframe_insertion
#. inframe_deletion
#. disruptive_inframe_deletion
#. 5_prime_UTR_truncation+exon_loss_variant
#. 3_prime_UTR_truncation+exon_loss
#. splice_region_variant
#. splice_branch_variant
#. stop_retained_variant
#. initiator_codon_variant
#. synonymous_variant
#. initiator_codon_variant+non_canonical_start_codon
#. stop_retained_variant
#. coding_sequence_variant
#. 5_prime_UTR_variant
#. 3_prime_UTR_variant
#. 5_prime_UTR_premature_start_codon_gain_variant
#. upstream_gene_variant
#. downstream_gene_variant
#. TF_binding_site_variant
#. regulatory_region_variant
#. miRNA
#. custom
#. sequence_feature
#. conserved_intron_variant
#. intron_variant
#. intragenic_variant
#. conserved_intergenic_variant
#. intergenic_region
#. coding_sequence_variant
#. non_coding_transcript_exon_variant
#. non_coding_transcript_variant
#. gene_variant
#. chromosome
