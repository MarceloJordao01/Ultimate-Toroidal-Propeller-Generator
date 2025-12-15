// toroidal_propeller.scad
include <BOSL2/std.scad>;

eps = 1/128;
$fn = 100;

// ============================================================
// Basics
// ============================================================
function lerp(a,b,u) = a + (b-a)*u;
function lerp2(p,q,u) = [lerp(p[0],q[0],u), lerp(p[1],q[1],u)];
function clamp(x,a,b) = (x<a)?a:((x>b)?b:x);

function rot2d_pt(p, ang) =
    let(c=cos(ang), s=sin(ang))
    [ p[0]*c - p[1]*s, p[0]*s + p[1]*c ];
function rot2d_path(path2d, ang) = [for(p=path2d) rot2d_pt(p, ang)];

function path2d_to3d(path2d) = [for(p=path2d) [p[0], p[1], 0]];

function vunit(v) = (norm(v) < 1e-9) ? [0,0,1] : v/norm(v);
function vcross(a,b) = cross(a,b);
function rot_about_axis(v, t, ang) =
    let(c=cos(ang), s=sin(ang))
    v*c + vcross(t,v)*s + t*(t*v)*(1-c);

function frame_matrix(p, n, b, t) =
    [[n[0],b[0],t[0],p[0]],
     [n[1],b[1],t[1],p[1]],
     [n[2],b[2],t[2],p[2]],
     [0,0,0,1]];

function path_tangent(path, i) =
    let(N=len(path),
        p0=(i==0)?path[0]:path[i-1],
        p1=(i==N-1)?path[N-1]:path[i+1])
    vunit(p1-p0);

function make_nb_from_tangent(t, up=[0,0,1]) =
    let(up2=(abs(t*up)>0.95)?[1,0,0]:up,
        b=vunit(vcross(t,up2)),
        n=vunit(vcross(b,t)))
    [n,b];

// ============================================================
// Keyframes (percent)
// ============================================================
function list_get_or_last(L, i, fallback=undef) =
    (L==undef || len(L)==0) ? fallback :
    (i < len(L)) ? L[i] : L[len(L)-1];

function kf_seg_index(pcts, u_pct) =
    let(n=len(pcts),
        hits=[for(j=[0:n-2]) if(u_pct>=pcts[j] && u_pct<=pcts[j+1]) j])
    (n<2)?0 : (len(hits)>0)?hits[len(hits)-1] : (u_pct<pcts[0])?0:(n-2);

function kf_local_u(pcts, j, u_pct) =
    let(a=pcts[j], b=pcts[j+1], d=b-a)
    (abs(d)<1e-9)?0:(u_pct-a)/d;

function kf_eval_seg(pcts, u_pct, nprofiles) =
    let(p=(pcts==undef||len(pcts)<2)?[0,100]:pcts,
        n=len(p),
        uu=clamp(u_pct,p[0],p[n-1]),
        j0=kf_seg_index(p,uu),
        j=clamp(j0,0,max(0,nprofiles-2)),
        v=kf_local_u(p,j0,uu))
    [j,v];

// angle wrap interpolation (evita pulo em -180/180)
function ang_norm(a) = let(x=(a+180)%360) x-180;
function ang_diff(a,b) = ang_norm(b-a);
function ang_lerp(a,b,u) = a + ang_diff(a,b)*u;

// ============================================================
// Path spline (Catmull–Rom 3D)
// ============================================================
function catmull_rom3d(u, P0, P1, P2, P3) =
    let(u2=u*u, u3=u2*u)
    [0.5*(2*P1[0] + (-P0[0]+P2[0])*u + (2*P0[0]-5*P1[0]+4*P2[0]-P3[0])*u2 + (-P0[0]+3*P1[0]-3*P2[0]+P3[0])*u3),
     0.5*(2*P1[1] + (-P0[1]+P2[1])*u + (2*P0[1]-5*P1[1]+4*P2[1]-P3[1])*u2 + (-P0[1]+3*P1[1]-3*P2[1]+P3[1])*u3),
     0.5*(2*P1[2] + (-P0[2]+P2[2])*u + (2*P0[2]-5*P1[2]+4*P2[2]-P3[2])*u2 + (-P0[2]+3*P1[2]-3*P2[2]+P3[2])*u3)];

function toroidal_path_spline_3d(
    t, hub_d, hub_height, blade_length, blade_offset,
    leadW_pct, trailW_pct, leadX_pct, trailX_pct,
    angle_A=60, angle_D=-240
) =
    let(
        leadW  = blade_length*(leadW_pct/100),
        trailW = blade_length*(trailW_pct/100),
        leadX  = blade_length*(leadX_pct/100),
        trailX = blade_length*(trailX_pct/100),
        r = hub_d*cos(30)/2,

        A = [ r*cos(60),  r*sin(60),   hub_height/2 + blade_offset/2 ],
        B = [ leadX,      +leadW,      hub_height/2 ],
        M = [ blade_length, 0,         hub_height/2 ],
        C = [ trailX,     -trailW,     hub_height/2 ],
        D = [ r*cos(60),  r*sin(-60),  hub_height/2 - blade_offset/2 ],

        scale_A = norm(B-A),
        scale_D = norm(D-C),

        Ta = [ scale_A*cos(angle_A), scale_A*sin(angle_A), 0 ],
        Td = [ scale_D*cos(angle_D), scale_D*sin(angle_D), 0 ],

        A0 = B - 2*Ta,
        D3 = C + 2*Td,

        segw=1/4,
        seg = (t<segw)?0 : (t<2*segw)?1 : (t<3*segw)?2 : 3,
        u   = (seg==0)?t/segw : (seg==1)?(t-segw)/segw : (seg==2)?(t-2*segw)/segw : (t-3*segw)/segw,

        P0=(seg==0)?A0 : (seg==1)?A : (seg==2)?B : M,
        P1=(seg==0)?A  : (seg==1)?B : (seg==2)?M : C,
        P2=(seg==0)?B  : (seg==1)?M : (seg==2)?C : D,
        P3=(seg==0)?M  : (seg==1)?C : (seg==2)?D : D3
    )
    catmull_rom3d(u,P0,P1,P2,P3);

function toroidal_path_points(
    steps, hub_d, hub_height, blade_length, blade_offset,
    leadW_pct, trailW_pct, leadX_pct, trailX_pct,
    t_start=0, t_end=1
) =
[
    for(i=[0:steps])
        let(u=i/steps, t=t_start + (t_end-t_start)*u)
        toroidal_path_spline_3d(t, hub_d, hub_height, blade_length, blade_offset,
                               leadW_pct, trailW_pct, leadX_pct, trailX_pct)
];

// ============================================================
// Profiles (NACA / flat / ellipse) + seam normalization
// ============================================================
function naca_digits_from_string(s) =
    [for(i=[0:len(s)-1]) ord(s[i]) - ord("0")];

function naca4_path(digits, chord=1, pts=40) =
    let(
        m=digits[0]/100,
        p=digits[1]/10,
        tt=(digits[2]*10+digits[3])/100,

        upper=[for(i=[0:pts]) let(
            xc=1-i/pts,
            yt=5*tt*(0.2969*sqrt(xc)-0.1260*xc-0.3516*xc*xc +0.2843*xc*xc*xc-0.1015*xc*xc*xc*xc),
            yc=(m==0)?0:(xc<p ? (m/(p*p))*(2*p*xc-xc*xc) : (m/((1-p)*(1-p)))*((1-2*p)+2*p*xc-xc*xc)),
            dyc=(m==0)?0:(xc<p ? (2*m/(p*p))*(p-xc) : (2*m/((1-p)*(1-p)))*(p-xc)),
            th=atan(dyc),
            xU=(xc-yt*sin(th))*chord,
            yU=(yc+yt*cos(th))*chord
        ) [xU,yU]],

        lower=[for(i=[0:pts]) let(
            xc=i/pts,
            yt=5*tt*(0.2969*sqrt(xc)-0.1260*xc-0.3516*xc*xc +0.2843*xc*xc*xc-0.1015*xc*xc*xc*xc),
            yc=(m==0)?0:(xc<p ? (m/(p*p))*(2*p*xc-xc*xc) : (m/((1-p)*(1-p)))*((1-2*p)+2*p*xc-xc*xc)),
            dyc=(m==0)?0:(xc<p ? (2*m/(p*p))*(p-xc) : (2*m/((1-p)*(1-p)))*(p-xc)),
            th=atan(dyc),
            xL=(xc+yt*sin(th))*chord,
            yL=(yc-yt*cos(th))*chord
        ) [xL,yL]],

        closed=concat(lower,upper)
    )
    [for(pp=closed) [(pp[0]-0.25*chord), pp[1]]];

function flat_plate_path(thickness=0.02, pts=40) =
    let(t=thickness/2)
    concat([for(i=[0:pts]) [(i/pts)-0.25, -t]],
           [for(i=[0:pts]) [(1-i/pts)-0.25,  t]]);

function ellipse_profile_path(thickness=0.12, pts=80) =
    let(t=thickness/2,
        upper=[for(i=[0:pts]) let(xc=1-i/pts, x1=2*xc-1, y=t*sqrt(max(0,1-x1*x1))) [xc-0.25, y]],
        lower=[for(i=[0:pts]) let(xc=i/pts,   x1=2*xc-1, y=-t*sqrt(max(0,1-x1*x1))) [xc-0.25, y]])
    concat(lower,upper);

function pad_or_trim(path, n) =
    (len(path)==n) ? path :
    (len(path)>n) ? [for(i=[0:n-1]) path[i]] :
    concat(path, [for(k=[0:(n-len(path)-1)]) path[len(path)-1]]);

function poly_area2(path) =
    sum([for(i=[0:len(path)-1])
        let(p=path[i], q=path[(i+1)%len(path)])
        (p[0]*q[1] - q[0]*p[1])]);

function reverse_path(path) = [for(i=[len(path)-1:-1:0]) path[i]];
function rotate_list(L,k) = (k<=0)?L:concat([for(i=[k:len(L)-1])L[i]],[for(i=[0:k-1])L[i]]);

function argmax_te_lower(path, tol=1e-6) =
    let(maxx=max([for(p=path)p[0]]),
        cand=[for(i=[0:len(path)-1]) if(abs(path[i][0]-maxx)<tol) i],
        besty=min([for(idx=cand) path[idx][1]]),
        best=[for(idx=cand) if(path[idx][1]==besty) idx])
    best[0];

function normalize_profile(path) =
    let(p0=path[0], pN=path[len(path)-1],
        base=(abs(p0[0]-pN[0])<1e-9 && abs(p0[1]-pN[1])<1e-9) ? [for(i=[0:len(path)-2])path[i]] : path,
        ccw=(poly_area2(base)<0) ? reverse_path(base) : base,
        k=argmax_te_lower(ccw))
    rotate_list(ccw,k);

function profile_path(profile_id, pts=80) =
    let(expected=2*(pts+1),
        raw =
            is_string(profile_id) ? naca4_path(naca_digits_from_string(profile_id), 1, pts) :
            (is_list(profile_id) && profile_id[0]=="flat")    ? flat_plate_path(profile_id[1], pts) :
            (is_list(profile_id) && profile_id[0]=="ellipse") ? ellipse_profile_path(profile_id[1], pts) :
            (is_list(profile_id) && profile_id[0]=="custom")  ? profile_id[1] :
            assert(false,"Perfil desconhecido / formato inválido."))
    pad_or_trim(normalize_profile(raw), expected);

// ============================================================
// Wing loft (skin)
// ============================================================
module toroidal_wing_skin_keyed(
    hub_d, hub_height, blade_length, blade_offset,
    leadW_pct, trailW_pct, leadX_pct, trailX_pct,
    profiles, profile_pcts,
    chords, attack_angles, chord_pivot_pcts,
    profile_pts=80, path_steps=64, path_portion=1.0,
    extra_twist_deg=-180, attack_zero_offset=90
){
    np = len(profiles);
    if (np < 2) echo("ERRO: profiles precisa ter >= 2.");
    else {
        path = toroidal_path_points(path_steps, hub_d, hub_height, blade_length, blade_offset,
                                    leadW_pct, trailW_pct, leadX_pct, trailX_pct, 0, path_portion);

        unit_profiles = [for(k=[0:np-1]) profile_path(profiles[k], profile_pts)];

        sections = [
            for(i=[0:len(path)-1])
                let(
                    u=(len(path)==1)?0:i/(len(path)-1),
                    u_pct=100*u,

                    seg=kf_eval_seg(profile_pcts, u_pct, np),
                    j=seg[0], v=seg[1],

                    chord0=list_get_or_last(chords, j, 6),
                    chord1=list_get_or_last(chords, j+1, chord0),
                    chord_u=lerp(chord0, chord1, v),

                    atk0=list_get_or_last(attack_angles, j, 0),
                    atk1=list_get_or_last(attack_angles, j+1, atk0),
                    atk_u=ang_lerp(atk0, atk1, v),

                    piv0=list_get_or_last(chord_pivot_pcts, j, 25)/100,
                    piv1=list_get_or_last(chord_pivot_pcts, j+1, piv0*100)/100,
                    piv_u=lerp(piv0, piv1, v),

                    p0=[for(p=unit_profiles[j])   [p[0]*chord_u, p[1]*chord_u]],
                    p1=[for(p=unit_profiles[j+1]) [p[0]*chord_u, p[1]*chord_u]],
                    prof_u=[for(k=[0:len(p0)-1]) lerp2(p0[k], p1[k], v)],

                    dx=(piv_u-0.25)*chord_u,
                    prof_piv=[for(p=prof_u) [p[0]+dx, p[1]]],

                    // (mantive sua convenção)
                    prof_rot=rot2d_path([for(p=prof_piv) [-p[0], p[1]]], atk_u + attack_zero_offset + 180),

                    t=path_tangent(path,i),
                    nb=make_nb_from_tangent(t,[0,0,1]),
                    n0=nb[0], b0=nb[1],

                    tw=u*extra_twist_deg,
                    n=rot_about_axis(n0,t,tw),
                    b=rot_about_axis(b0,t,tw),

                    T=frame_matrix(path[i], n, b, t)
                )
                apply(T, path2d_to3d(prof_rot))
        ];

        skin(sections, slices=1, caps=true);
    }
}

// ============================================================
// Main module
// ============================================================
module toroidal_propeller(
    blades=3,
    blade_length=68,
    blade_offset=-6,

    hub_height=6,
    hub_d=16,
    hub_screw_d=5.5,
    hub_notch_height=0,
    hub_notch_d=0,

    leading_edge_blade_width=18,
    trailing_edge_blade_width=18,
    leading_edge_blade_xoffset=25,
    trailing_edge_blade_xoffset=25,

    profiles=["2412",["ellipse",0.5],"2412"],
    profile_pcts=[0,50,100],

    chords=[6,2,6],
    attack_angles=[15,0,-10],
    chord_pivot_pcts=[0,0,0],

    profile_pts=80,
    path_steps=64,
    path_portion=1.0,
    extra_twist_deg=-180,
    attack_zero_offset=90
){
    difference() {
        union() {
            for(a=[0:blades-1])
                rotate([0,0,a*360/blades])
                    toroidal_wing_skin_keyed(
                        hub_d, hub_height, blade_length, blade_offset,
                        leading_edge_blade_width, trailing_edge_blade_width,
                        leading_edge_blade_xoffset, trailing_edge_blade_xoffset,
                        profiles, profile_pcts,
                        chords, attack_angles, chord_pivot_pcts,
                        profile_pts, path_steps, path_portion,
                        extra_twist_deg, attack_zero_offset
                    );

            rotate([0,0,30]) cylinder(d=hub_d, h=hub_height, $fn=6);
        }

        translate([0,0,-eps]) cylinder(d=hub_screw_d, h=hub_height+2*eps);

        if (hub_notch_height>0 && hub_notch_d>0)
            translate([0,0,-eps]) cylinder(d=hub_notch_d, h=hub_notch_height+eps);
    }
}

// default
toroidal_propeller();
