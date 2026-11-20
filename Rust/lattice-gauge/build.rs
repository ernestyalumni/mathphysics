// Builds cuda/su2_update.cu with nvcc only when the `cuda` feature is enabled. The default build
// stays std-only and needs no CUDA toolkit.
fn main() {
    println!("cargo:rerun-if-changed=cuda/su2_update.cu");
    println!("cargo:rerun-if-changed=build.rs");
    if std::env::var_os("CARGO_FEATURE_CUDA").is_none() {
        return;
    }
    let out = std::env::var("OUT_DIR").unwrap();
    let cuda = std::env::var("CUDA_HOME").unwrap_or_else(|_| "/usr/local/cuda".into());
    let arch = std::env::var("LG_CUDA_ARCH").unwrap_or_else(|_| "sm_86".into()); // RTX 3060
    let obj = format!("{out}/su2_update.o");
    let status = std::process::Command::new(format!("{cuda}/bin/nvcc"))
        .args(["-O3", "-fmad=false", "-Xcompiler", "-fPIC", &format!("-arch={arch}"), "-c", "cuda/su2_update.cu", "-o", &obj])
        .status()
        .expect("failed to run nvcc");
    assert!(status.success(), "nvcc failed");
    let status = std::process::Command::new("ar").args(["crs", &format!("{out}/liblg_cuda.a"), &obj]).status().expect("ar");
    assert!(status.success(), "ar failed");
    println!("cargo:rustc-link-search=native={out}");
    println!("cargo:rustc-link-lib=static=lg_cuda");
    println!("cargo:rustc-link-search=native={cuda}/lib64");
    println!("cargo:rustc-link-lib=dylib=cudart");
    println!("cargo:rustc-link-lib=dylib=stdc++");
}
