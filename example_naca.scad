// naca_segment_linear_yz_rot.scad
//
// Segmento de asa NACA indo de p0 até p1,
// com perfil inicial e final diferentes,
// perfis no plano YZ e rotacionados em Z.
//
//  - p0 = [1,0,0]
//  - p1 = [5,1,0]

include <src/BOSL2/std.scad>;
include <src/naca_profilles.scad>;

//------------------------------------------------------------
// Helpers
//------------------------------------------------------------

// Interpolação linear
function lerp(a, b, t) = a + (b - a) * t;

// Rotaciona um ponto 3D em torno do eixo Z
function rotate_z(pt, ang) =
    let(
        c = cos(ang),
        s = sin(ang)
    )
    [
        c*pt[0] - s*pt[1],  // x'
        s*pt[0] + c*pt[1],  // y'
        pt[2]               // z'
    ];

// Converte perfil 2D (x,y) → 3D no plano YZ, rotaciona em Z e translada
function translate_profile_3d_yz_rot(profile2d, pos, rot_z_deg=30) =
    [
        for (p = profile2d)
        let(
            // (x, y) do perfil → (0, x, y) no espaço 3D
            pt0 = [0, p[0], p[1]],
            // rotaciona em torno do eixo Z
            pt1 = rotate_z(pt0, rot_z_deg)
        )
        [ pt1[0] + pos[0], pt1[1] + pos[1], pt1[2] + pos[2] ]
    ];

// Interpola ponto a ponto entre dois perfis 2D (mesmo número de pontos)
function interp_profile_2d(root2d, tip2d, t) =
    [
        for (i = [0 : len(root2d) - 1])
        [
            lerp(root2d[i][0], tip2d[i][0], t),
            lerp(root2d[i][1], tip2d[i][1], t)
        ]
    ];

// Path linear de p0 até p1
function line_path(t, p0=[1,0,0], p1=[5,1,0]) =
    [
        lerp(p0[0], p1[0], t),
        lerp(p0[1], p1[1], t),
        lerp(p0[2], p1[2], t)
    ];

//------------------------------------------------------------
// Módulo principal: segmento NACA entre dois pontos
//------------------------------------------------------------

module naca_segment_linear(
    p0 = [1, 0, 0],
    p1 = [5, 1, 0],

    root_code    = "2412",
    tip_code     = "0012",

    chord_root   = 1.0,
    chord_tip    = 0.5,

    n_pts        = 160,   // pontos ao longo da corda
    steps_path   = 20,    // seções ao longo do caminho

    rot_z_deg    = 30     // rotação em torno de Z para todos os perfis
) {

    // Perfis 2D de início e fim (em XY, saída padrão da NACA)
    root2d = naca4_profile(root_code, chord=chord_root, n=n_pts);
    tip2d  = naca4_profile(tip_code,  chord=chord_tip,  n=n_pts);

    // Gera lista de perfis 3D ao longo do caminho linear
    profiles3d = [
        for (i = [0 : steps_path])
        let(
            t     = i / steps_path,
            sec2d = interp_profile_2d(root2d, tip2d, t),
            pos   = line_path(t, p0, p1)
        )
        translate_profile_3d_yz_rot(sec2d, pos, rot_z_deg)
    ];

    // Loft com BOSL2
    skin(
        profiles = profiles3d,
        slices   = 0,
        caps     = true
    );
}

//------------------------------------------------------------
// Chamada de exemplo
//------------------------------------------------------------

naca_segment_linear(
    p0         = [1, 0, 0],
    p1         = [5, 1, 0],
    root_code  = "2412",
    tip_code   = "0012",
    chord_root = 1.0,
    chord_tip  = 0.5,
    n_pts      = 160,
    steps_path = 20,
    rot_z_deg  = 30
);
