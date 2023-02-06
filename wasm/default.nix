{ stdenv, fetchFromGitHub, callPackage, lib }:

let
  csound7Sources = fetchFromGitHub {
    owner = "csound";
    repo = "csound";
    rev = "c9b5c8e939cc47156dc412a11a28d7f0d9632a61";
    sha256 = "sha256-kWUeoXOS7fedYQCIaz/Kcb/vk2rJbCzryApQdHe67Lo=";
  };

  csound6Sources = fetchFromGitHub {
    owner = "csound";
    repo = "csound";
    rev = "74c91c04bdf20386f52aae13da36494a6b5ca4be";
    sha256 = "sha256-jYQVoBMBDrWptPRVzgKnj7r3iNPp6A06fW0j8Gk7n+U=";
  };

  wasi-sdk = callPackage "${csound6Sources}/wasm/src/wasi-sdk.nix" { };
  csound-wasm = callPackage "${csound6Sources}/wasm/src/csound.nix" {};

in stdenv.mkDerivation {
  name = "hypercurve-wasm-plugin";
  buildInputs = [ csound-wasm ];

  src =  ../.;
  dontStrip = true;

  buildPhase = ''

    echo "Compile hypercurve.wasm"
    ${wasi-sdk}/bin/clang++  \
      -fPIC -std=c++17 -Wno-deprecated-declarations \
      -D__wasi__=1 \
      --target=wasm32-unknown-emscripten \
      --sysroot=${wasi-sdk}/share/wasi-sysroot \
      -DUSE_DOUBLE=1 \
      -I${wasi-sdk}/share/wasi-sysroot/include \
      -I${csound-wasm}/include \
      -I../src \
      -I../src/asciiplot \
      -I../src/fpng/src \
      -c  \
      csound_opcode/csound_hypercurve.cpp \
      src/asciiplot/asciiplotter.cpp \
      src/fpng/src/fpng.cpp

    ${wasi-sdk}/bin/wasm-ld \
      -z stack-size=128 \
      -pie --no-entry \
      --experimental-pic \
      --import-table \
      --import-memory \
      --export=__wasm_call_ctors \
      --export=csoundModuleInit \
      --export-if-defined=csoundModuleCreate \
      --export-if-defined=csoundModuleDestroy \
      -L${wasi-sdk}/share/wasi-sysroot/lib/wasm32-unknown-emscripten \
      -L${csound-wasm}/lib -lcsound-wasm \
      -lc -lc++ -lc++abi \
      *.o -o hypercurve.wasm

 # -lrt -lutil -lxnet -lresolv -lc-printscan-long-double \
 #      -lwasi-emulated-getpid -lwasi-emulated-signal -lwasi-emulated-mman -lwasi-emulated-process-clocks \
  '';

  installPhase = ''
    mkdir -p $out/lib
    cp ./hypercurve.wasm $out/lib
  '';
}
