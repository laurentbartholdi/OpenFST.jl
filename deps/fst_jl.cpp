#include "jlcxx/jlcxx.hpp"
#include "jlcxx/tuple.hpp"

#include <fst/fstlib.h>

static std::stringstream errormsg;

fst::StdVectorFst testerror(fst::StdVectorFst &f)
{
  if (f.Properties(fst::kError, true)) {
    std::string s = errormsg.str();
    std::string prefix = "LEVEL(FST_FLAGS_fst_error_fatal ? base_logging::FATAL : base_logging::ERROR):";
    int pos;
    while ((pos = s.find(prefix)) != std::string::npos)
      s.replace(pos, prefix.length(), "(OpenFST)");
    while (!s.empty() && s.back() == '\n')
      s.pop_back();
    jl_error(s.c_str());
  }
    
  return f;
}

// standard visitor class to return DFS ordering.
class OrderVisitor {
 public:
  using Arc = typename fst::StdArc;
  using StateId = typename Arc::StateId;

  // If acyclic, order[i] gives the topological position of StateId i;
  // otherwise it is unchanged. acyclic_ will be true iff the FST has no
  // cycles. The caller retains ownership of the state order vector.
  OrderVisitor(std::vector<StateId> *order) : order_(order) {}

  void InitVisit(const fst::Fst<Arc> &fst) { visit_.clear(); }

  bool InitState(StateId s, StateId) { visit_.push_back(s); return true; }

  constexpr bool TreeArc(StateId, const Arc &) const { return true; }

  bool BackArc(StateId, const Arc &) { return true; }

  constexpr bool ForwardOrCrossArc(StateId, const Arc &) const { return true; }

  void FinishState(StateId s, StateId, const Arc *) const { }

  void FinishVisit() {
    order_->clear();
    for (StateId s = 0; s < visit_.size(); ++s) {
      order_->push_back(fst::kNoStateId);
    }
    for (StateId s = 0; s < visit_.size(); ++s) {
      (*order_)[visit_[s]] = s;
    }
  }

 private:
  std::vector<StateId> *order_, visit_;
};
  
/* put a transducer in canonical form: first sort the arcs, and then renumber the states with start=0 and visited along DFS.
 */
fst::StdVectorFst canonize(fst::StdVectorFst &f) {
  fst::ArcSort(&f, fst::StdILabelCompare());
  std::vector<fst::StdArc::StateId> order;
  OrderVisitor visitor(&order);
  fst::DfsVisit(f, &visitor);
  fst::StateSort(&f, order);
  
  return f;
}

// visitor class to recursively delete leaves
class LeafVisitor {
 public:
  using Arc = typename fst::StdArc;
  using StateId = typename Arc::StateId;

  LeafVisitor(std::vector<StateId> *del) : delete_(del) {}

  void InitVisit(const fst::Fst<Arc> &f) { leaf.clear(); }

  bool InitState(StateId s, StateId) { leaf.resize(std::max(leaf.size(),(size_t)s+1)); leaf[s] = true; return true; }

  constexpr bool TreeArc(StateId, const Arc &) const { return true; }

  bool BackArc(StateId s, const Arc &) { leaf[s] = false; return true; }

  bool ForwardOrCrossArc(StateId s, const Arc &a) { leaf[s] &= leaf[a.nextstate]; return true; }

  void FinishState(StateId s, StateId parent, const Arc *) { if (parent != fst::kNoStateId) leaf[parent] &= leaf[s]; }

  void FinishVisit() {
    delete_->clear();
    for (StateId s = 0; s < leaf.size(); ++s)
      if (leaf[s])
	delete_->push_back(s);
  }

 private:
  std::vector<StateId> *delete_;
  std::vector<char> leaf; // bool would make it read-only in my compiler
};
  
/* remove from a transducer all states that don't accept an infinite word.
 */
fst::StdVectorFst omegawords(fst::StdVectorFst &f) {
  fst::Connect(&f);
  std::vector<fst::StdArc::StateId> del;
  LeafVisitor visitor(&del);
  fst::DfsVisit(f, &visitor);
  f.DeleteStates(del);
  
  return f;
}

// nn2n(x,y) = (x+y)*(x+y+1)÷2+x
// n2nn(z) = (s = (Int(floor(sqrt(8z+1)))-1)÷2; (z-s*(s+1)÷2,s*(s+3)÷2-z))

constexpr unsigned cantor(unsigned x, unsigned y) { return (x+y)*(x+y+1)/2 + x; }
std::pair<unsigned,unsigned> rotnac(unsigned z) {
  unsigned s = (sqrt(8*z+1.1)-1)/2;
  return std::pair<unsigned,unsigned>(z-s*(s+1)/2,s*(s+3)/2-z);
}

struct CantorMapper {
  using Arc = typename fst::StdArc;
  fst::MapFinalAction FinalAction() const { return fst::MAP_NO_SUPERFINAL; }
  fst::MapSymbolsAction InputSymbolsAction() const { return fst::MAP_NOOP_SYMBOLS; }
  fst::MapSymbolsAction OutputSymbolsAction() const { return fst::MAP_NOOP_SYMBOLS; }
  unsigned long Properties(unsigned long props) const { return (props & ~(fst::kAcceptor | fst::kNotAcceptor)) | fst::kAcceptor; }
  Arc operator()(const Arc &arc) { unsigned z = cantor(arc.ilabel,arc.olabel); return Arc(z,z,arc.weight,arc.nextstate); }
};

struct RotnacMapper {
  using Arc = typename fst::StdArc;
  constexpr explicit RotnacMapper(void) : acceptor(true) { }
  fst::MapFinalAction FinalAction() const { return fst::MAP_NO_SUPERFINAL; }
  fst::MapSymbolsAction InputSymbolsAction() const { return fst::MAP_NOOP_SYMBOLS; }
  fst::MapSymbolsAction OutputSymbolsAction() const { return fst::MAP_NOOP_SYMBOLS; }
  unsigned long Properties(unsigned long props) const { return (props & ~(fst::kAcceptor | fst::kNotAcceptor)) | (acceptor ? fst::kAcceptor : fst::kNotAcceptor); }
  Arc operator()(const Arc &arc) { auto z = rotnac(arc.ilabel); acceptor &= (z.first == z.second); return Arc(z.first,z.second,arc.weight,arc.nextstate); }
private:
  bool acceptor;
};

/* a variant of "encode": replace a transducer by an acceptor */
fst::StdVectorFst transducer2acceptor(fst::StdVectorFst &f) {
  CantorMapper mapper;
  fst::ArcMap(&f, mapper);
  return f;
}

fst::StdVectorFst transducer2acceptor_st(fst::StdVectorFst &f) {
  auto t = fst::SymbolTable();
  for (auto &a : *f.InputSymbols())
    for (auto &b : *f.OutputSymbols())
      t.AddSymbol(a.Symbol()+":"+b.Symbol(), cantor(a.Label(),b.Label()));
  f.SetInputSymbols(&t);
  f.SetOutputSymbols(&t);
  return f;
}

/* a variant of "decode": replace an acceptor by a transducer */
fst::StdVectorFst acceptor2transducer(fst::StdVectorFst &f) {
  RotnacMapper mapper;
  fst::ArcMap(&f, mapper);
  //  f.SetProperties(fst::kAcceptor | fst::kNotAcceptor, fst::kNotAcceptor);
  return f;
}

fst::StdVectorFst acceptor2transducer_st(fst::StdVectorFst &f) {
  auto t = fst::SymbolTable(), u = fst::SymbolTable();
  for (auto &a : *f.InputSymbols())
    for (auto &b : *f.OutputSymbols()) {
      auto colon = a.Symbol().find(":");
      auto z = rotnac(a.Label());
      t.AddSymbol(a.Symbol().substr(0,colon), z.first);
      u.AddSymbol(a.Symbol().substr(colon+1), z.second);
    }
  f.SetInputSymbols(&t);
  f.SetOutputSymbols(&u);
  return f;
}

fst::StdVectorFst state(fst::StdVectorFst &f, int in, int out) {
  fst::StdVectorFst g(f);

  if (f.Start() != fst::kNoStateId)
    for (fst::ArcIterator<fst::StdVectorFst> aiter(f, f.Start()); !aiter.Done(); aiter.Next()) {
      auto a = aiter.Value();
      if (a.ilabel == in && a.olabel == out) {
	g.SetStart(a.nextstate);
	return g;
      }
    }
  g.DeleteStates();
  return g;
}
  
JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  FST_FLAGS_fst_error_fatal = false;
  FST_FLAGS_outputstream = &errormsg;
  
  mod.add_type<fst::SymbolTable>("SymbolTable")
    .method("addtable!", &fst::SymbolTable::AddTable)
    .method("availablekey", &fst::SymbolTable::AvailableKey)
    .method("checksum", &fst::SymbolTable::LabeledCheckSum)
    .method("getnthkey", &fst::SymbolTable::GetNthKey)
    .method("numsymbols", &fst::SymbolTable::NumSymbols)
    .method("removesymbol", &fst::SymbolTable::RemoveSymbol);
  mod.method("addsymbol!", [](fst::SymbolTable &s, std::string name) { return s.AddSymbol(name); });
  mod.method("addsymbol!", [](fst::SymbolTable &s, std::string name, int id) { return s.AddSymbol(name, id); });
  mod.method("find", [](fst::SymbolTable &s, std::string name) { return s.Find(name); });
  mod.method("_find", [](fst::SymbolTable &s, int id) { return s.Find(id); });
  mod.method("member", [](fst::SymbolTable &s, std::string name) { return s.Member(name); });
  mod.method("member", [](fst::SymbolTable &s, int id) { return s.Member(id); });
  
  mod.add_type<fst::TropicalWeight>("TropicalWeight");

  mod.add_type<fst::StdArc>("StdArc")
    .constructor<int, int, float, int>()
    .method("splat", [](fst::StdArc &a) { return std::make_tuple(a.ilabel, a.olabel, a.weight.Value(), a.nextstate); });
  
  mod.add_type<fst::StdVectorFst>("StdVectorFst")
    .constructor<fst::StdVectorFst>()
    .method("properties", &fst::StdVectorFst::Properties)
    .method("_setinputsymbols!", &fst::StdVectorFst::SetInputSymbols)
    .method("_inputsymbols", &fst::StdVectorFst::InputSymbols)
    .method("_setoutputsymbols!", &fst::StdVectorFst::SetOutputSymbols)
    .method("_outputsymbols", &fst::StdVectorFst::OutputSymbols)
    .method("numstates", &fst::StdVectorFst::NumStates)
    .method("addstate!", &fst::StdVectorFst::AddState)
    .method("setstart!", &fst::StdVectorFst::SetStart)
    .method("_setfinal!", &fst::StdVectorFst::SetFinal)
    .method("start",  &fst::StdVectorFst::Start)
    .method("_final",  &fst::StdVectorFst::Final)
    ;
  mod.method("addarc!", [](fst::StdVectorFst &f, int from, fst::StdArc arc) { f.AddArc(from, arc); });
  mod.method("_addarc!", [](fst::StdVectorFst &f, int from, int to, int ilabel, int olabel, float weight) { f.AddArc(from, fst::StdArc(ilabel, olabel, weight, to)); });
  mod.method("setfinal!", [](fst::StdVectorFst &f, int to, float weight) { f.SetFinal(to, weight); });
  mod.method("final", [](fst::StdVectorFst &f, int from) { return f.Final(from).Value(); });
  mod.method("isexpanded", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kNotAcceptor, true); });
  mod.method("ismutable", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kMutable, true); });
  mod.method("iserror", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kError, true); });
  mod.method("isacceptor", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kAcceptor, true); });
  mod.method("isideterministic", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kIDeterministic, true); });
  mod.method("isodeterministic", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kODeterministic, true); });
  mod.method("hasepsilons", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kEpsilons, true); });
  mod.method("hasiepsilons", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kIEpsilons, true); });
  mod.method("hasoepsilons", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kOEpsilons, true); });
  mod.method("isilabelsorted", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kILabelSorted, true); });
  mod.method("isolabelsorted", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kOLabelSorted, true); });
  mod.method("isweighted", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kWeighted, true); });
  mod.method("iscyclic", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kCyclic, true); });
  mod.method("isinitialcyclic", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kInitialCyclic, true); });
  mod.method("istopsorted", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kTopSorted, true); });
  mod.method("isaccessible", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kAccessible, true); });
  mod.method("iscoaccessible", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kCoAccessible, true); });
  mod.method("isstring", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kString, true); });
  mod.method("hasweightedcycles", [](fst::StdVectorFst &f) -> bool { return f.Properties(fst::kWeightedCycles, true); });

  mod.method("testerror", testerror);

  mod.method("_read", [](std::vector<unsigned char> &s) { std::stringstream ss(std::string(s.begin(),s.end())); return fst::StdVectorFst::Read(ss, fst::FstReadOptions()); });
  mod.method("_write", [](fst::StdVectorFst &f) { std::stringstream ss; f.Write(ss, fst::FstWriteOptions()); std::string s = ss.str(); return std::vector<unsigned char>(s.c_str(),s.c_str()+s.length()); });
  mod.method("_read", [](std::string s) { return fst::StdVectorFst::Read(s); });
  mod.method("_write", [](std::string s, fst::StdVectorFst &f) { return f.Write(s); });
  
  mod.add_type<fst::StateIterator<fst::StdVectorFst>>("StateIterator")
    .constructor<fst::StdVectorFst>()
    .method("Done", &fst::StateIterator<fst::StdVectorFst>::Done)
    .method("Value", &fst::StateIterator<fst::StdVectorFst>::Value)
    .method("Next", &fst::StateIterator<fst::StdVectorFst>::Next);
  
  mod.add_type<fst::ArcIterator<fst::StdVectorFst>>("ArcIterator")
    .constructor<fst::StdVectorFst,int>()
    .method("Done", &fst::ArcIterator<fst::StdVectorFst>::Done)
    .method("Value", &fst::ArcIterator<fst::StdVectorFst>::Value)
    .method("Next", &fst::ArcIterator<fst::StdVectorFst>::Next);

  mod.add_type<fst::EncodeMapper<fst::StdArc>>("EncodeMapper");
  
  mod.method("_arcsort_i!", [](fst::StdVectorFst &f) { fst::ArcSort(&f, fst::StdILabelCompare()); return testerror(f); });
  mod.method("_arcsort_o!", [](fst::StdVectorFst &f) { fst::ArcSort(&f, fst::StdOLabelCompare()); return testerror(f); });

  mod.method("_closure_s!", [](fst::StdVectorFst &f) { fst::Closure(&f, fst::CLOSURE_STAR); return testerror(f); });
  mod.method("_closure_p!", [](fst::StdVectorFst &f) { fst::Closure(&f, fst::CLOSURE_PLUS); return testerror(f); });

  mod.method("_compose", [](fst::StdVectorFst &f, fst::StdVectorFst &g) { fst::StdVectorFst h; fst::Compose(f, g, &h); return testerror(h); });

  mod.method("concat!",  [](fst::StdVectorFst &f, fst::StdVectorFst &g) { fst::Concat(f, &g); return testerror(g); });

  mod.method("connect!",  [](fst::StdVectorFst &f) { fst::Connect(&f); return f; });

  mod.method("determinize", [](fst::StdVectorFst &f) { fst::StdVectorFst g; fst::Determinize(f, &g); return testerror(g); });

  mod.method("difference", [](fst::StdVectorFst &f, fst::StdVectorFst &g) { fst::StdVectorFst h; fst::Difference(f, g, &h); return testerror(h); });
  
  mod.method("disambiguate", [](fst::StdVectorFst &f) { fst::StdVectorFst g; fst::Disambiguate(f, &g); return testerror(g); });

  mod.method("_epsnormalize_i", [](fst::StdVectorFst &f) { fst::StdVectorFst g; fst::EpsNormalize(f, &g, fst::EPS_NORM_INPUT); return testerror(g); });
  mod.method("_epsnormalize_o", [](fst::StdVectorFst &f) { fst::StdVectorFst g; fst::EpsNormalize(f, &g, fst::EPS_NORM_OUTPUT); return testerror(g); });

  mod.method("equivalent", [](fst::StdVectorFst &f, fst::StdVectorFst &g) { return fst::Equivalent(f, g); });

  mod.method("_intersect", [](fst::StdVectorFst &f, fst::StdVectorFst &g) { fst::StdVectorFst h; fst::Intersect(f, g, &h); return testerror(h); });

  mod.method("inv!", [](fst::StdVectorFst &f) { fst::Invert(&f); return testerror(f); });
  
  mod.method("isomorphic", [](fst::StdVectorFst &f, fst::StdVectorFst &g) { return fst::Isomorphic(f, g); });

  mod.method("minimize!", [](fst::StdVectorFst &f, bool allow_nondet) { fst::Minimize(&f, (fst::StdMutableFst *) nullptr, fst::kDelta, allow_nondet); return testerror(f); });
  mod.method("minimize_rt!", [](fst::StdVectorFst &f, bool allow_nondet) { fst::StdVectorFst g; fst::Minimize(&f, &g, fst::kDelta, allow_nondet); testerror(f); return std::make_tuple(f,g); });

  mod.method("_project_i", [](fst::StdVectorFst &f) { fst::StdVectorFst g; fst::Project(f, &g, fst::ProjectType::INPUT); return testerror(g); });
  mod.method("_project_o", [](fst::StdVectorFst &f) { fst::StdVectorFst g; fst::Project(f, &g, fst::ProjectType::OUTPUT); return testerror(g); });

  mod.method("_push_i", [](fst::StdVectorFst &f) { fst::StdVectorFst g; fst::Push<fst::StdArc, fst::REWEIGHT_TO_INITIAL>(f, &g, fst::kPushLabels); return testerror(g); });
  mod.method("_push_f", [](fst::StdVectorFst &f) { fst::StdVectorFst g; fst::Push<fst::StdArc, fst::REWEIGHT_TO_FINAL>(f, &g, fst::kPushLabels); return testerror(g); });
 
  mod.method("rmepsilon!", [](fst::StdVectorFst &f) { fst::RmEpsilon(&f); return testerror(f); });

  mod.method("synchronize", [](fst::StdVectorFst &f) { fst::StdVectorFst g; fst::Synchronize(f, &g); return testerror(g); });

  mod.method("topsort!", [](fst::StdVectorFst &f) { fst::TopSort(&f); return testerror(f); });

  mod.method("_union!", [](fst::StdVectorFst &f, fst::StdVectorFst &g) { fst::Union(&f, g); return testerror(f); });

  mod.method("verify", [](fst::StdVectorFst &f) { return fst::Verify(f); });

  mod.method("encode!", [](fst::StdVectorFst &f) { fst::EncodeMapper<fst::StdArc> encoder(fst::kEncodeLabels, fst::ENCODE); fst::Encode(&f, &encoder); return encoder; });
  mod.method("encode!", [](fst::StdVectorFst &f, fst::EncodeMapper<fst::StdArc>  &encoder) { fst::Encode(&f, &encoder); return f; });
  mod.method("decode!", [](fst::StdVectorFst &f, fst::EncodeMapper<fst::StdArc>  &encoder) { fst::Decode(&f, encoder); return f; });

  mod.method("canonize!", &canonize);
  mod.method("omegawords!", &omegawords);
  mod.method("state", &state);
  mod.method("transducer2acceptor!", &transducer2acceptor);
  mod.method("transducer2acceptor_st!", &transducer2acceptor_st);
  mod.method("acceptor2transducer!", &acceptor2transducer);
  mod.method("acceptor2transducer_st!", &acceptor2transducer_st);
  
  // Base methods
  mod.set_override_module(jl_base_module);
  mod.method("copy", [](fst::SymbolTable &s) { return *s.Copy(); });
  mod.method("copy!", [](fst::SymbolTable &s, fst::SymbolTable &t) { return s = t; });

  mod.method("copy", [](fst::StdVectorFst &f) { return *f.Copy(); });
  mod.method("copy!", [](fst::StdVectorFst &f, fst::StdVectorFst &g) { return f = g; });
  mod.method("==", [](fst::StdVectorFst &f, fst::StdVectorFst &g) { return fst::Equal(f, g); });
  mod.method("reverse", [](fst::StdVectorFst &f) { fst::StdVectorFst g; fst::Reverse(f, &g); return testerror(g); });
  mod.unset_override_module();
}
