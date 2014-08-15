helpers do
  def get_gene(gene, organism = "Hsa/jun2011")
    gene = Translation.job(:translate, nil, :genes => [gene], :organism => organism).exec.first
    Gene.setup(gene, "Ensembl Gene ID", organism)
  end

  def get_genes(genes, organism = "Hsa/jun2011")
    gene = Translation.job(:translate, nil, :genes => genes, :organism => organism).exec
    Gene.setup(gene, "Ensembl Gene ID", organism)
  end

  def get_mutation(mutation, organism = "Hsa/jun2011")
    GenomicMutation.setup(mutation, nil, organism, true)
  end
end
