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


helpers do

  def load_user_studies
    if $user_studies.nil?

      public_studies = Study.studies.select do |study| study.users.nil? or study.users.include? 'public' end
      $user_studies = Hash.new{ public_studies.dup }
      $user_studies.extend IndiferentHash

      groups = Rbbt.etc.web_user_groups.exists? ? (Rbbt.etc.web_user_groups.yaml || {}) : {}
      Study.studies.each do |study|
        next if study.metadata.nil? or study.metadata[:hide]
        users = study.users.nil? ? [:public] : study.users.collect{|user| groups.include?(user)? groups[user] : user}.flatten.uniq
        users.each do |user|
          $user_studies[user.to_sym] += [study]
        end
      end
    end
    $user_studies
  end

  def user_studies
    load_user_studies if $user_studies.nil?
    $user_studies 
  end
end

before do
  parts = request.path_info.split "/"
  parts.shift
  case
  when parts[0..1] == %w(entity Study)
    study = parts[2]
  when parts[0..1] == %w(entity_action Study)
    study = parts[3]
  else
    study = nil
  end

  authorize! if study

  #if study and ( not user_studies.include? user or not user_studies[user].include? study )
  if study and not user_studies[user].include? study 
    halt 401, "We are sorry, but it seems that you do not have credentials to explore this study: #{ study }"
  end
end

