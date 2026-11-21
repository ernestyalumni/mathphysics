import Lean

/-!
# Statement export for the statement checker

`#ym_export foo` prints the elaborated type of the constant `foo` as one line of JSON, prefixed by
`YM_EXPORT `. Binders are instantiated with free variables, so the Rust side never sees de Bruijn
indices. Every binder records its user name, whether its type is a proposition, and its binder info.
Our own definitions (not Lean's or Mathlib's) that the statement mentions are exported with their values;
structures in that namespace are exported with their constructor type, so a hypothesis such as
`Spectrum lam` can be read as the conjunction of its fields.
-/

open Lean Meta Elab Command

namespace YMAudit

private def binderInfoString : BinderInfo → String
  | .default => "default"
  | .implicit => "implicit"
  | .strictImplicit => "strict"
  | .instImplicit => "inst"

/-- JSON for an expression, with binders opened as free variables. -/
partial def toJson (e : Expr) : MetaM Json := do
  match e with
  | .bvar i => return Json.mkObj [("k", "bvar"), ("i", i)]
  | .fvar id => do
      let name := (← id.getDecl).userName
      return Json.mkObj [("k", "fvar"), ("id", toString id.name), ("name", toString name)]
  | .mvar _ => return Json.mkObj [("k", "mvar")]
  | .sort l => return Json.mkObj [("k", "sort"), ("prop", l.isZero)]
  | .const n _ => return Json.mkObj [("k", "const"), ("n", toString n)]
  | .app .. => do
      let args ← e.getAppArgs.mapM toJson
      return Json.mkObj [("k", "app"), ("f", ← toJson e.getAppFn), ("args", Json.arr args)]
  | .lam n t b bi => binder "lam" n t b bi
  | .forallE n t b bi => do
      let selfProp ← (try isProp e catch _ => pure false)
      let j ← binder "pi" n t b bi
      return j.setObjVal! "self_prop" (Json.bool selfProp)
  | .letE n t v b _ => do
      let tj ← toJson t
      let vj ← toJson v
      withLetDecl n t v fun x => do
        return Json.mkObj [("k", "let"), ("name", toString n), ("id", toString x.fvarId!.name),
          ("type", tj), ("value", vj), ("body", ← toJson (b.instantiate1 x))]
  | .lit (.natVal v) => return Json.mkObj [("k", "nat"), ("v", toString v)]
  | .lit (.strVal s) => return Json.mkObj [("k", "str"), ("v", s)]
  | .mdata _ e => toJson e
  | .proj s i e => return Json.mkObj [("k", "proj"), ("s", toString s), ("i", i), ("e", ← toJson e)]
where
  binder (kind : String) (n : Name) (t b : Expr) (bi : BinderInfo) : MetaM Json := do
    let tj ← toJson t
    let isPropType ← (try isProp t catch _ => pure false)
    withLocalDecl n bi t fun x => do
      let body := b.instantiate1 x
      return Json.mkObj [("k", kind), ("name", toString n), ("id", toString x.fvarId!.name),
        ("bi", binderInfoString bi), ("prop", isPropType), ("dep", b.hasLooseBVars),
        ("type", tj), ("body", ← toJson body)]

/-- Is `n` defined in this file or in one of our own modules (not Lean, Std or Mathlib)? -/
def isOurs (env : Environment) (n : Name) : Bool :=
  match env.getModuleIdxFor? n with
  | none => true
  | some idx =>
      let m := env.header.moduleNames[idx.toNat]!
      !((`Mathlib).isPrefixOf m || (`Init).isPrefixOf m || (`Lean).isPrefixOf m || (`Std).isPrefixOf m
        || (`Batteries).isPrefixOf m || (`Aesop).isPrefixOf m || (`Qq).isPrefixOf m || (`Plausible).isPrefixOf m)

/-- Our own constants mentioned by `e`. -/
def localConstants (env : Environment) (e : Expr) : Array Name :=
  e.getUsedConstants.filter (isOurs env)

/-- Export `root` and, transitively, the `YMFormal` definitions and structures it uses. -/
def exportConstant (root : Name) : MetaM Json := do
  let info ← getConstInfo root
  let env ← getEnv
  let mut todo := localConstants env info.type
  let mut seen : NameSet := {}
  let mut defs : Array Json := #[]
  while h : todo.size > 0 do
    let n := todo.back
    todo := todo.pop
    if seen.contains n then continue
    seen := seen.insert n
    let ci ← getConstInfo n
    match ci with
    | .defnInfo d =>
        defs := defs.push (Json.mkObj [("name", toString n), ("kind", "def"), ("value", ← toJson d.value)])
        todo := todo ++ localConstants env d.value
    | .inductInfo d =>
        if isStructure (← getEnv) n then
          let ctor ← getConstInfo d.ctors.head!
          defs := defs.push (Json.mkObj [("name", toString n), ("kind", "structure"),
            ("params", d.numParams), ("ctor", ← toJson ctor.type)])
          todo := todo ++ localConstants env ctor.type
    | _ => pure ()
  return Json.mkObj [("name", toString root), ("type", ← toJson info.type), ("defs", Json.arr defs)]

syntax (name := ymExport) "#ym_export " ident : command

@[command_elab ymExport] def elabYmExport : CommandElab := fun stx => do
  let id := stx[1]
  let names ← liftCoreM <| realizeGlobalConstWithInfos id
  for n in names do
    let json ← liftTermElabM <| exportConstant n
    logInfo m!"YM_EXPORT {json.compress}"

end YMAudit
