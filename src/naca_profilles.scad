// src/naca_profilles.scad
// Geração de perfis NACA 4 dígitos.
//
// Saída: lista de [x,y] em volta do aerofólio,
// começando no bordo de fuga inferior,
// indo até o bordo de ataque,
// e voltando pelo extradorso até o bordo de fuga superior.
//
// Uso típico:
//
//     include <src/naca_profilles.scad>;
//
//     pts = naca4_profile("2412", chord=1, n=160);
//     polygon(pts);
//


//-------------------------
// Funções auxiliares simples
//-------------------------

// Converte um caractere numérico em inteiro (0–9)
function _digit(str, idx) = ord(str[idx]) - ord("0");

// Lê dois dígitos consecutivos como inteiro (00–99)
function _two_digits(str, idx) =
    10 * _digit(str, idx) + _digit(str, idx+1);


//-------------------------
// Função principal NACA 4 dígitos
//-------------------------
//
// code  : string de 4 dígitos, ex: "2412", "0012"
// chord : corda em unidades do modelo
// n     : número de subdivisões ao longo da corda
//
// Retorna: lista de [x,y]
//
function naca4_profile(code="2412", chord=1, n=80) =
    let(
        // Parâmetros NACA MPXX
        m = _digit(code, 0) / 100,       // camber
        p = _digit(code, 1) / 10,        // posição da camber
        t = _two_digits(code, 2) / 100,  // espessura relativa

        // Distribuição LINEAR em x pra evitar “deformação esquisita”
        xs = [ for (i = [0:n]) chord * i/n ]
    )
    concat(
        // Intradorso (baixo): TE -> LE
        [
            for (i = [n : -1 : 0])
            let(
                x  = xs[i],
                xc = x/chord,

                yt = 5*t * (
                        0.2969*sqrt(xc)
                      - 0.1260*xc
                      - 0.3516*xc*xc
                      + 0.2843*xc*xc*xc
                      - 0.1015*xc*xc*xc*xc
                    ),

                yc = (m == 0 || p == 0) ? 0 :
                     ( xc <= p
                       ?  m/(p*p)            * (2*p*xc - xc*xc)
                       :  m/((1-p)*(1-p))    * ((1 - 2*p) + 2*p*xc - xc*xc)
                     ),

                dyc_dx = (m == 0 || p == 0) ? 0 :
                         ( xc <= p
                           ?  2*m/(p*p)         * (p - xc)
                           :  2*m/((1-p)*(1-p)) * (p - xc)
                         ),

                // OpenSCAD usa graus em atan/sin/cos
                th = atan(dyc_dx)
            )
            [ x + yt*sin(th), yc - yt*cos(th) ]
        ],

        // Extradorso (cima): LE -> TE
        [
            for (i = [0 : n])
            let(
                x  = xs[i],
                xc = x/chord,

                yt = 5*t * (
                        0.2969*sqrt(xc)
                      - 0.1260*xc
                      - 0.3516*xc*xc
                      + 0.2843*xc*xc*xc
                      - 0.1015*xc*xc*xc*xc
                    ),

                yc = (m == 0 || p == 0) ? 0 :
                     ( xc <= p
                       ?  m/(p*p)            * (2*p*xc - xc*xc)
                       :  m/((1-p)*(1-p))    * ((1 - 2*p) + 2*p*xc - xc*xc)
                     ),

                dyc_dx = (m == 0 || p == 0) ? 0 :
                         ( xc <= p
                           ?  2*m/(p*p)         * (p - xc)
                           :  2*m/((1-p)*(1-p)) * (p - xc)
                         ),

                th = atan(dyc_dx)
            )
            [ x - yt*sin(th), yc + yt*cos(th) ]
        ]
    );


//-------------------------
// (Opcional) Módulo de teste 2D
//-------------------------
module demo_naca2d() {
    pts = naca4_profile("2412", chord=1, n=160);
    polygon(pts);
}
