module Study
  def self.watson(study)
    study_info(study)[:watson]
  end

  def has_genotypes?
    @has_genotypes ||= samples.select{|s| s.has_genotype? }.any?
  end

  attr_accessor :watson
  def watson
    @watson  = study_info[:watson] if @watson.nil?
    @watson
  end

  def cohort
    @cohort ||= samples.collect do |sample| 
      next unless sample.has_genotype?
      genomic_mutations = sample.genomic_mutations
      GenomicMutation.setup(genomic_mutations, sample, organism, watson)
    end.compact
  end

  def self.value_stream(stream)
    _s = TSV.traverse stream, :into => :stream, :type => :flat do |key,elems|
      next if elems.nil? or elems.empty?
      elems * "\n"
    end
    CMD.cmd('sort -u', :in => _s, :pipe => true)
  end

  property :metagenotype => :single do
    GenomicMutation.setup(cohort.flatten.sort, self + " metagenotype", organism, watson)
  end

  property :get_genes => :single do |type=nil|
    @sample_gene_matches ||= subset(:sample_genes, :source => :all, :target => :all)
    case type
    when :recurrent, "recurrent", "Recurrent"
      recurrent = Misc.counts(@sample_gene_matches.select_by(:info){|info| info["affected"] == "true"}.target).select{|g,c| c > 1 }.collect{|g,c| g }
      @sample_gene_matches.target_entity.uniq.subset(recurrent)
    when nil
      @sample_gene_matches.target_entity.uniq
    else
      @sample_gene_matches.select_by(:info){|info| info[type.to_s] == "true" }.target_entity.uniq
    end
  end
end
