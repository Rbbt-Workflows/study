helpers do
  def get_gene(gene, organism = Organism.default_code("Hsa"))
    gene = Translation.job(:translate, nil, :genes => [gene], :organism => organism).exec.first
    Gene.setup(gene, "Ensembl Gene ID", organism)
  end

  def get_genes(genes, organism = Organism.default_code("Hsa"))
    gene = Translation.job(:translate, nil, :genes => genes, :organism => organism).exec
    Gene.setup(gene, "Ensembl Gene ID", organism)
  end

  def get_mutation(mutation, organism = Organism.default_code("Hsa"))
    GenomicMutation.setup(mutation, nil, organism, true)
  end
end

Misc.bootstrap(Study.studies, 3) do |study|
  begin
    Study.setup(study).knowledge_base
  rescue
    Log.warn $!.message
  end
end

TSV.traverse Study.studies, :bar => "Registering study KBs" do |study|
  begin
    KnowledgeBaseRESTHelpers.add_syndication study, Study.setup(study).knowledge_base
  rescue
    Log.warn $!.message
  end
end

