module Study

  def self.mutation_classification(mis, type)
    if mis.any?
      case
      when mis.select{|mi| mi =~ /Frame/}.any?
        if type == "INS"
          return "Frame_Shift_Ins"
        else
          return "Frame_Shift_Del"
        end
      when mis.select{|mi| mi =~ /ENSP.*:([A-Z])\d+([A-Z])/ && $1 != $2 }.any?
        return "Missense_Mutation"
      when mis.select{|mi| mi =~ /ENSP.*:[A-Z]\d+\*/}.any?
        return "Nonsense_Mutation"
      when mis.select{|mi| mi =~ /ENSP.*:\*\d+[A-Z]/ }.any?
        return "Nonstop_Mutation"
      when mis.select{|mi| mi =~ /ENSP.*:Indel/}.any?
        if type == "INS"
          return "In_Frame_Ins"
        else
          return "In_Frame_Del"
        end
      when mis.select{|mi| mi =~ /ENST.*:UTR3/ }.any?
        return "3'UTR"
      when mis.select{|mi| mi =~ /ENST.*:UTR5/ }.any?
        return "5'UTR"
      when mis.select{|mi| mi =~ /ENSP.*:([A-Z\*])\d+([A-Z\*])/ && $1 == $2 }.any?
        return "Silent"
      else
        raise "Unkown: #{mis * ", "}"
      end
    else
      classification = "IGR"
    end
  end

  MAF_FIELDS = %w(
    Entrez_Gene_Id
    Center
    NCBI_Build
    Chromosome
    Start_Position
    End_Position
    Strand
    Variant_Classification
    Variant_Type
    Reference_Allele
    Tumor_Seq_Allele1
    Tumor_Seq_Allele2
    dbSNP_RS
    dbSNP_Val_Status
    Tumor_Sample_Barcode
    Matched_Norm_Sample_Barc
    Match_Norm_Seq_Allele1
    Match_Norm_Seq_Allele2
    Tumor_Validation_Allele1
    Tumor_Validation_Allele2
    Match_Norm_Validation_Allele1
    Match_Norm_Validation_Allele2
    Verification_Status
    Validation_Status
    Mutation_Status
    Sequencing_Phase
    Sequence_Source
    Validation_Method
    Score
    BAM_File
    Sequencer
    Tumor_Sample_UUID
    Matched_Norm_Sample_UUID
  )


  MISSING = "MISSING"
  dep :genomic_mutations
  dep :mutation_incidence
  dep :genomic_mutation_consequence
  dep Sequence, :reference, :positions => :genomic_mutations, :organism => :organism
  task :maf_file => :text do

    ensg2name = Organism.identifiers(organism).index :target => "Associated Gene Name", :fields => ["Ensembl Gene ID"], :persist => true, :unnamed => true
    ensp2ensg = Organism.transcripts(organism).index :target => "Ensembl Gene ID", :fields => ["Ensembl Protein ID"], :persist => true, :unnamed => true
    enst2ensg = Organism.transcripts(organism).index :target => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"], :persist => true, :unnamed => true
    gene_strand = Organism.gene_positions(organism).tsv :fields => ["Strand"], :type => :single, :persist => true, :unnamed => true
    pasted = TSV.paste_streams([step(:reference), step(:mutation_incidence), step(:genomic_mutation_consequence)])
    dumper = TSV::Dumper.new :key_field => "Hugo_Symbol", :fields => MAF_FIELDS, :type => :list
    dumper.init
    TSV.traverse pasted, :type => :double, :into => dumper do |mutation,values|
      reference, samples, mis = values

      mutation = mutation.first if Array === mutation
      reference = reference.first if Array === reference

      mutation_parts = mutation.split(":")

      result = []
      result.extend MultipleResult

      protein = mis.any? ? mis.first.split(":").first : nil
      entrez = ensp2ensg[protein] || enst2ensg[protein]
      gene = ensg2name[entrez] || "Unknown"

      chr = mutation_parts[0]
      start = mutation_parts[1]
      allele = mutation_parts[2]

      type = allele.length == 1 ? "SNP" : (allele[0] == '+' ? "INS" : "DEL")

      classification = Study.mutation_classification(mis, type)

      center = ""
      build = ""

      eend = start.to_i + allele.length - 1
      strand = gene_strand[entrez] == 1 ? '+' : '-'


      reference = reference
      allele2 = allele
      rs = ""
      rs_validation = ""
      norm_barcode = ""
      norm_allele = ""
      norm_allele2 = ""
      validation_allele = ""
      validation_allele2 = ""
      norm_validation_allele = ""
      norm_validation_allele2 = ""
      verification = ""
      validation = "Untested"
      status = "Somatic"
      phase = ""
      source = "WGS"
      validation_method = "none"
      score = ""
      bam = ""
      sequencer = "UNKOWN"
      norm_uuid = ""

      samples.each do |sample|
        barcode = sample
        uuid = sample

        values = []
        values << entrez
        values << center
        values << build
        values << chr
        values << start
        values << eend
        values << strand
        values << classification
        values << type
        values << reference
        values << allele
        values << allele2
        values << rs
        values << rs_validation
        values << barcode
        values << norm_barcode
        values << norm_allele
        values << norm_allele2
        values << validation_allele
        values << validation_allele2
        values << norm_validation_allele
        values << norm_validation_allele2
        values << verification
        values << validation
        values << status
        values << phase
        values << source
        values << validation_method
        values << score
        values << bam
        values << sequencer
        values << uuid
        values << norm_uuid

        result << [gene, values]
      end
      result
    end
  end

  #Workflow.require_workflow "MutSig"
  #if defined? MutSig
  #  dep :maf_file
  #  dep MutSig, :analysis, :maf_file => :maf_file
  #  task :mut_sig => :tsv do
  #    TSV.get_stream step(:analysis)
  #  end
  #end
end
