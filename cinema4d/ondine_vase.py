# -*- coding: utf-8 -*-
"""
PICSEL  |  ONDINE VASE  |  Cinema 4D Python (Script Manager)
Cible : Cinema 4D 2024 uniquement (Python 3.9). Unites : centimetres.

Vase ceramique contemporain, silhouette en S legerement asymetrique :
  - ventre decale d'un cote, col qui penche de l'autre (equilibre dynamique),
  - ouverture legerement elliptique,
  - trois fines gorges sur l'epaule (rythme par repetition),
  - exterieur craie mat, interieur emaille vert profond (contraste),
  - pied biseaute, levre arrondie.

Le maillage est genere par anneaux (100 % quads, sauf 2 fans de triangles
au centre du fond), pas de subdivision necessaire.
Lancer : Script Manager > Nouveau > coller > Executer (Ctrl+Entree).
"""

import math
import c4d

# =============================================================================
# PARAMETRES (modifiables)
# =============================================================================
ROOT_NAME = "PICSEL_Ondine_Vase"

HEIGHT = 32.0            # hauteur totale (cm)
WALL = 0.5               # epaisseur de paroi (cm)
FLOOR_THICKNESS = 1.0    # epaisseur du fond (cm)
FOOT_BEVEL = 0.5         # rayon du biseau du pied (cm)
SEGMENTS = 96            # divisions radiales (multiple de 4 conseille)

# Profil exterieur : (hauteur y, rayon r) - interpole par Catmull-Rom
PROFILE_CTRL = [
    (FOOT_BEVEL, 5.4),
    (2.5, 7.4),
    (7.0, 9.8),
    (12.5, 10.5),
    (18.0, 8.8),
    (23.0, 5.6),
    (27.0, 3.9),
    (HEIGHT - WALL * 0.5, 4.3),
]

# Gorges sur l'epaule : (hauteur, profondeur, largeur sigma)
GROOVES = [(20.0, 0.14, 0.32), (21.4, 0.14, 0.32), (22.8, 0.14, 0.32)]

# Asymetrie controlee
NECK_LEAN_X = 2.4        # decalage du col en +X (cm)
NECK_LEAN_Z = 0.0
BELLY_SHIFT_X = -0.9     # decalage du ventre en -X (cm), cree le "S"
MOUTH_ELLIPSE = 0.10     # 0 = rond ; 0.10 = ellipse douce a l'ouverture

# Resolution du maillage le long du profil
STEP_COARSE = 0.9
STEP_FINE = 0.28
RIM_STEPS = 10           # pair : le point de la levre est au sommet
FOOT_STEPS = 5

# Sens des normales. Si le vase apparait "a l'envers" (Phong sombre, faces
# retournees en mode backface culling), passer a True.
FLIP_NORMALS = False

# Materiaux
COL_EXTERIOR = (0.90, 0.88, 0.84)     # craie chaude
COL_INTERIOR = (0.07, 0.17, 0.13)     # vert bouteille profond
COL_FLOOR = (0.80, 0.78, 0.74)
ENABLE_MICRO_BUMP = True              # leger grain de ceramique (cosmetique)

# Rendu Etsy (portrait 4:5)
RENDER_W, RENDER_H = 2400, 3000


# =============================================================================
# OUTILS
# =============================================================================
def log(msg):
    print("[Ondine] " + msg)


def const(name):
    """Retourne la constante c4d.<name> ou None (avec avertissement)."""
    value = getattr(c4d, name, None)
    if value is None:
        log("ATTENTION : constante c4d.%s absente dans cette version (ignoree)." % name)
    return value


def set_param(node, name, value):
    """Affecte node[c4d.<name>] = value ; avertit sans planter si indisponible."""
    pid = const(name)
    if pid is None:
        return False
    try:
        node[pid] = value
    except (TypeError, AttributeError, KeyError) as exc:
        log("ATTENTION : parametre %s non affecte (%s)." % (name, exc))
        return False
    return True


def clamp01(x):
    return max(0.0, min(1.0, x))


def smoothstep(x):
    x = clamp01(x)
    return x * x * (3.0 - 2.0 * x)


# =============================================================================
# GEOMETRIE (Python pur, sans API c4d)
# =============================================================================
def catmull_rom_curve(ctrl, per_segment=48):
    """Polyligne dense passant par les points de controle (uniform Catmull-Rom)."""
    pts = list(ctrl)
    # points fantomes aux extremites (reflexion) pour definir les tangentes
    first = (2 * pts[0][0] - pts[1][0], 2 * pts[0][1] - pts[1][1])
    last = (2 * pts[-1][0] - pts[-2][0], 2 * pts[-1][1] - pts[-2][1])
    pts = [first] + pts + [last]
    out = []
    for i in range(1, len(pts) - 2):
        p0, p1, p2, p3 = pts[i - 1], pts[i], pts[i + 1], pts[i + 2]
        for s in range(per_segment):
            t = s / float(per_segment)
            t2, t3 = t * t, t * t * t
            pt = []
            for k in (0, 1):
                pt.append(0.5 * ((2 * p1[k]) + (-p0[k] + p2[k]) * t
                                 + (2 * p0[k] - 5 * p1[k] + 4 * p2[k] - p3[k]) * t2
                                 + (-p0[k] + 3 * p1[k] - 3 * p2[k] + p3[k]) * t3))
            out.append((pt[0], pt[1]))
    out.append(ctrl[-1])
    return out


def resample(poly, step_fn):
    """Re-echantillonne une polyligne (y, r) a pas variable (pas = step_fn(y))."""
    lengths = [0.0]
    for i in range(1, len(poly)):
        lengths.append(lengths[-1] + math.hypot(poly[i][0] - poly[i - 1][0],
                                                poly[i][1] - poly[i - 1][1]))
    total = lengths[-1]
    out = [poly[0]]
    s = step_fn(poly[0][0])
    idx = 1
    while s < total - 1e-6:
        while lengths[idx] < s:
            idx += 1
        seg = lengths[idx] - lengths[idx - 1]
        u = (s - lengths[idx - 1]) / seg if seg > 0 else 0.0
        y = poly[idx - 1][0] + u * (poly[idx][0] - poly[idx - 1][0])
        r = poly[idx - 1][1] + u * (poly[idx][1] - poly[idx - 1][1])
        out.append((y, r))
        s += step_fn(y)
    if math.hypot(out[-1][0] - poly[-1][0], out[-1][1] - poly[-1][1]) < 0.15 and len(out) > 1:
        out[-1] = poly[-1]
    else:
        out.append(poly[-1])
    return out


def groove_offset(y):
    return sum(d * math.exp(-((y - yc) / w) ** 2) for yc, d, w in GROOVES)


def fine_step(y):
    near = any(abs(y - yc) < 1.1 for yc, _, _ in GROOVES)
    return STEP_FINE if near else STEP_COARSE


def build_profile():
    """
    Profil (r, y) du bas-centre vers l'exterieur, le haut, la levre, l'interieur,
    puis le centre du fond interieur. Retourne (points, index_milieu_levre).
    Le sens de parcours garantit des normales sortantes.
    """
    base = catmull_rom_curve(PROFILE_CTRL)

    outer = resample(base, fine_step)
    outer = [(y, r - groove_offset(y)) for (y, r) in outer]
    inner = resample(base, lambda y: STEP_COARSE)

    prof = [(0.0, 0.0)]                                   # centre du dessous

    # Pied biseaute : quart de cercle de (rb-b, 0) a (rb, b)
    rb = PROFILE_CTRL[0][1]
    b = FOOT_BEVEL
    for i in range(FOOT_STEPS + 1):
        a = math.radians(-90.0 + 90.0 * i / FOOT_STEPS)
        prof.append((rb - b + b * math.cos(a), b + b * math.sin(a)))

    # Paroi exterieure (le 1er point = fin du biseau, deja present)
    for (y, r) in outer[1:]:
        prof.append((r, y))

    # Levre : demi-cercle de l'exterieur vers l'interieur
    rm = outer[-1][1]
    y_top = outer[-1][0]
    rad = WALL * 0.5
    cx, cy = rm - rad, y_top
    mid_index = None
    for i in range(1, RIM_STEPS + 1):
        a = math.pi * i / RIM_STEPS
        prof.append((cx + rad * math.cos(a), cy + rad * math.sin(a)))
        if i == RIM_STEPS // 2:
            mid_index = len(prof) - 1

    # Paroi interieure (descente), decalee de WALL vers l'axe
    y_min = FLOOR_THICKNESS + 0.8
    for (y, r) in reversed(inner[:-1]):
        if y >= y_min:
            prof.append((r - WALL, y))

    # Fond interieur
    r_last = prof[-1][0]
    prof.append((max(r_last - 0.5, 0.3), FLOOR_THICKNESS))
    prof.append((0.0, FLOOR_THICKNESS))
    return prof, mid_index


def deform(r, y, theta):
    """Position 3D d'un point du profil a l'angle theta, avec asymetrie."""
    w = smoothstep((y - 6.0) / (HEIGHT - 6.0))            # 0 en bas, 1 en haut
    belly = math.exp(-((y - 11.0) / 6.0) ** 2)
    cx = NECK_LEAN_X * w * w + BELLY_SHIFT_X * belly
    cz = NECK_LEAN_Z * w * w
    sx = 1.0 + MOUTH_ELLIPSE * w
    sz = 1.0 - MOUTH_ELLIPSE * w
    return (r * math.cos(theta) * sx + cx, y, r * math.sin(theta) * sz + cz)


def build_mesh_data():
    """
    Retourne (points, polygons, uvs, interior_polys)
      points   : liste de (x, y, z)
      polygons : liste de tuples d'indices (3 ou 4)
      uvs      : liste (par polygone) de tuples (u, v) de meme longueur
      interior_polys : indices des polygones de l'emaille interieur
    """
    prof, mid_index = build_profile()
    seg = SEGMENTS

    # v = abscisse curviligne normalisee du profil
    acc = [0.0]
    for i in range(1, len(prof)):
        acc.append(acc[-1] + math.hypot(prof[i][0] - prof[i - 1][0],
                                        prof[i][1] - prof[i - 1][1]))
    total = acc[-1]

    points = []
    ring_start = []
    for (r, y) in prof:
        ring_start.append(len(points))
        if r <= 1e-9:
            points.append(deform(0.0, y, 0.0))
        else:
            for j in range(seg):
                points.append(deform(r, y, 2.0 * math.pi * j / seg))

    def pid(k, j):
        if prof[k][0] <= 1e-9:
            return ring_start[k]
        return ring_start[k] + (j % seg)

    polygons, uvs, interior = [], [], []
    for k in range(len(prof) - 1):
        v0, v1 = acc[k] / total, acc[k + 1] / total
        pole0 = prof[k][0] <= 1e-9
        pole1 = prof[k + 1][0] <= 1e-9
        for j in range(seg):
            u0, u1 = j / float(seg), (j + 1) / float(seg)
            if pole0:
                poly = (pid(k, j), pid(k + 1, j), pid(k + 1, j + 1))
                uv = ((u0, v0), (u0, v1), (u1, v1))
            elif pole1:
                poly = (pid(k, j), pid(k + 1, j), pid(k, j + 1))
                uv = ((u0, v0), (u0, v1), (u1, v0))
            else:
                poly = (pid(k, j), pid(k + 1, j), pid(k + 1, j + 1), pid(k, j + 1))
                uv = ((u0, v0), (u0, v1), (u1, v1), (u1, v0))
            if FLIP_NORMALS:
                poly = tuple(reversed(poly))
                uv = tuple(reversed(uv))
            if k >= mid_index:
                interior.append(len(polygons))
            polygons.append(poly)
            uvs.append(uv)
    return points, polygons, uvs, interior


# =============================================================================
# OBJETS CINEMA 4D
# =============================================================================
def make_null(name, parent=None, doc=None):
    null = c4d.BaseObject(c4d.Onull)
    null.SetName(name)
    if doc is not None:
        doc.InsertObject(null, parent)
    return null


def make_vase_object():
    points, polygons, uvs, interior = build_mesh_data()
    obj = c4d.PolygonObject(len(points), len(polygons))
    obj.SetName("Ondine_Vase_Body")
    obj.SetAllPoints([c4d.Vector(*p) for p in points])
    for i, poly in enumerate(polygons):
        if len(poly) == 4:
            obj.SetPolygon(i, c4d.CPolygon(poly[0], poly[1], poly[2], poly[3]))
        else:
            obj.SetPolygon(i, c4d.CPolygon(poly[0], poly[1], poly[2]))

    # UVW cylindriques (pour grain, decals, exports vers d'autres logiciels)
    uvw = c4d.UVWTag(len(polygons))
    for i, uv in enumerate(uvs):
        vec = [c4d.Vector(u, v, 0.0) for (u, v) in uv]
        if len(vec) == 3:
            vec.append(vec[2])
        uvw.SetSlow(i, vec[0], vec[1], vec[2], vec[3])
    obj.InsertTag(uvw)

    # Selection de polygones pour l'emaille interieur
    sel_tag = c4d.BaseTag(c4d.Tpolygonselection)
    sel_tag.SetName("SEL_Interior")
    bs = sel_tag.GetBaseSelect()
    for i in interior:
        bs.Select(i)
    obj.InsertTag(sel_tag)

    # Lissage Phong (angle limite pour garder pied et levre nets)
    phong = c4d.BaseTag(c4d.Tphong)
    set_param(phong, "PHONGTAG_PHONG_ANGLELIMIT", True)
    set_param(phong, "PHONGTAG_PHONG_ANGLE", math.radians(50.0))
    obj.InsertTag(phong)

    obj.Message(c4d.MSG_UPDATE)
    return obj


# =============================================================================
# MATERIAUX
# =============================================================================
def configure_reflection(mat, roughness, strength):
    """Couche de reflexion GGX dielectrique (API reflectance R20+)."""
    if not hasattr(mat, "AddReflectionLayer"):
        log("ATTENTION : AddReflectionLayer indisponible, reflexion non configuree.")
        return
    set_param(mat, "MATERIAL_USE_REFLECTION", True)
    # Retire les couches par defaut pour partir d'une base propre
    for _ in range(16):
        if mat.GetReflectionLayerCount() == 0:
            break
        mat.RemoveReflectionLayerIndex(0)
    layer = mat.AddReflectionLayer()
    base = layer.GetDataID()
    for name, value in (
        ("REFLECTION_LAYER_MAIN_DISTRIBUTION", getattr(c4d, "REFLECTION_DISTRIBUTION_GGX", None)),
        ("REFLECTION_LAYER_MAIN_VALUE_ROUGHNESS", roughness),
        ("REFLECTION_LAYER_MAIN_VALUE_REFLECTION", strength),
        ("REFLECTION_LAYER_FRESNEL_MODE", getattr(c4d, "REFLECTION_FRESNEL_DIELECTRIC", None)),
    ):
        pid = getattr(c4d, name, None)
        if pid is None or value is None:
            log("ATTENTION : %s indisponible, valeur par defaut conservee." % name)
            continue
        mat[base + pid] = value


def add_micro_bump(mat, strength):
    if not hasattr(c4d, "Xnoise"):
        log("ATTENTION : shader Noise indisponible, bump ignore.")
        return
    shader = c4d.BaseShader(c4d.Xnoise)
    set_param(shader, "SLA_NOISE_SCALE", 0.04)
    mat.InsertShader(shader)
    set_param(mat, "MATERIAL_USE_BUMP", True)
    set_param(mat, "MATERIAL_BUMP_SHADER", shader)
    set_param(mat, "MATERIAL_BUMP_STRENGTH", strength)


def make_material(name, color, roughness, strength, bump=0.0):
    mat = c4d.BaseMaterial(c4d.Mmaterial)
    mat.SetName(name)
    set_param(mat, "MATERIAL_USE_COLOR", True)
    set_param(mat, "MATERIAL_COLOR_COLOR", c4d.Vector(*color))
    set_param(mat, "MATERIAL_USE_SPECULAR", False)
    configure_reflection(mat, roughness, strength)
    if bump > 0.0 and ENABLE_MICRO_BUMP:
        add_micro_bump(mat, bump)
    mat.Update(True, True)
    return mat


def make_backdrop_material(name, color):
    """Materiau lumineux uniforme pour le fond (pas d'ombrage parasite)."""
    mat = c4d.BaseMaterial(c4d.Mmaterial)
    mat.SetName(name)
    set_param(mat, "MATERIAL_USE_COLOR", False)
    set_param(mat, "MATERIAL_USE_LUMINANCE", True)
    set_param(mat, "MATERIAL_LUMINANCE_COLOR", c4d.Vector(*color))
    set_param(mat, "MATERIAL_USE_SPECULAR", False)
    set_param(mat, "MATERIAL_USE_REFLECTION", False)
    mat.Update(True, True)
    return mat


def assign_material(obj, mat, restriction=None):
    tag = c4d.BaseTag(c4d.Ttexture)
    tag[c4d.TEXTURETAG_MATERIAL] = mat
    if restriction:
        tag[c4d.TEXTURETAG_RESTRICTION] = restriction
    obj.InsertTag(tag)
    return tag


# =============================================================================
# PRESENTATION : sol, fond, lumieres, camera
# =============================================================================
def look_at_matrix(pos, target, up=None):
    up = up or c4d.Vector(0, 1, 0)
    fwd = (target - pos).GetNormalized()
    right = up.Cross(fwd).GetNormalized()
    true_up = fwd.Cross(right)
    return c4d.Matrix(pos, right, true_up, fwd)


def make_floor():
    floor = c4d.BaseObject(c4d.Oplane)
    floor.SetName("Studio_Floor")
    set_param(floor, "PRIM_PLANE_WIDTH", 600.0)
    set_param(floor, "PRIM_PLANE_HEIGHT", 600.0)
    set_param(floor, "PRIM_PLANE_SUBW", 1)
    set_param(floor, "PRIM_PLANE_SUBH", 1)
    set_param(floor, "PRIM_AXIS", getattr(c4d, "PRIM_AXIS_YP", None))
    return floor


def make_area_light(name, pos, target, intensity, size, color=(1.0, 0.97, 0.93)):
    light = c4d.BaseObject(c4d.Olight)
    light.SetName(name)
    set_param(light, "LIGHT_TYPE", getattr(c4d, "LIGHT_TYPE_AREA", None))
    set_param(light, "LIGHT_COLOR", c4d.Vector(*color))
    set_param(light, "LIGHT_BRIGHTNESS", intensity)
    set_param(light, "LIGHT_AREADETAILS_SIZEX", size[0])
    set_param(light, "LIGHT_AREADETAILS_SIZEY", size[1])
    set_param(light, "LIGHT_SHADOWTYPE", getattr(c4d, "LIGHT_SHADOWTYPE_AREA", None))
    light.SetMg(look_at_matrix(pos, target))
    return light


def make_camera(name, pos, target, focal_mm):
    cam = c4d.BaseObject(c4d.Ocamera)
    cam.SetName(name)
    set_param(cam, "CAMERA_FOCUS", focal_mm)
    cam.SetMg(look_at_matrix(pos, target))
    return cam


# =============================================================================
# ASSEMBLAGE
# =============================================================================
def remove_previous(doc):
    """Rend le script rejouable : supprime l'ancienne version."""
    old = doc.SearchObject(ROOT_NAME)
    while old is not None:
        old.Remove()
        old = doc.SearchObject(ROOT_NAME)
    for name in ("MAT_Ceramic_Chalk_Exterior", "MAT_Glaze_Deep_Green_Interior",
                 "MAT_Studio_Floor", "MAT_Studio_Backdrop"):
        mat = doc.SearchMaterial(name)
        while mat is not None:
            mat.Remove()
            mat = doc.SearchMaterial(name)


def build_scene(doc):
    remove_previous(doc)
    doc.StartUndo()

    root = make_null(ROOT_NAME, None, doc)
    model = make_null("MODEL", root, doc)
    stage = make_null("STAGE", root, doc)
    lights = make_null("LIGHTS", root, doc)
    camnull = make_null("CAMERA", root, doc)

    # --- Modele
    vase = make_vase_object()
    doc.InsertObject(vase, model)

    # --- Materiaux (dans le Material Manager, prefixe MAT_)
    mat_ext = make_material("MAT_Ceramic_Chalk_Exterior", COL_EXTERIOR, 0.42, 0.35, bump=0.03)
    mat_int = make_material("MAT_Glaze_Deep_Green_Interior", COL_INTERIOR, 0.12, 0.60)
    mat_floor = make_material("MAT_Studio_Floor", COL_FLOOR, 0.55, 0.15)
    mat_back = make_backdrop_material("MAT_Studio_Backdrop", COL_FLOOR)
    for m in (mat_ext, mat_int, mat_floor, mat_back):
        doc.InsertMaterial(m)

    # Exterieur d'abord, interieur ensuite (le tag restreint prend le dessus)
    assign_material(vase, mat_ext)
    assign_material(vase, mat_int, restriction="SEL_Interior")

    # --- Sol + fond
    floor = make_floor()
    doc.InsertObject(floor, stage)
    assign_material(floor, mat_floor)
    background = c4d.BaseObject(c4d.Obackground)
    background.SetName("Studio_Backdrop")
    doc.InsertObject(background, stage)
    assign_material(background, mat_back)

    # --- Lumieres (cible : milieu du vase)
    target = c4d.Vector(0.5, 14.0, 0.0)
    key = make_area_light("Light_Key_SoftBox", c4d.Vector(-95, 120, -85), target, 1.0, (90, 90))
    fill = make_area_light("Light_Fill", c4d.Vector(110, 55, -60), target, 0.35, (120, 60),
                           color=(0.92, 0.95, 1.0))
    rim = make_area_light("Light_Rim", c4d.Vector(30, 70, 110), target, 0.6, (60, 100))
    for light in (key, fill, rim):
        doc.InsertObject(light, lights)

    # --- Camera 85 mm, legerement en plongee
    cam = make_camera("Cam_Hero_85mm", c4d.Vector(18, 38, -118), c4d.Vector(0.8, 15.5, 0.0), 85.0)
    doc.InsertObject(cam, camnull)

    # --- Rendu portrait 4:5 (format Etsy)
    rd = doc.GetActiveRenderData()
    if rd is not None:
        set_param(rd, "RDATA_XRES", float(RENDER_W))
        set_param(rd, "RDATA_YRES", float(RENDER_H))

    # Camera de scene
    bd = doc.GetActiveBaseDraw()
    if bd is not None:
        bd.SetSceneCamera(cam)

    doc.AddUndo(c4d.UNDOTYPE_NEW, root)
    doc.EndUndo()
    return vase


def main():
    doc = c4d.documents.GetActiveDocument()
    if doc is None:
        log("ERREUR : aucun document actif.")
        return
    vase = build_scene(doc)
    c4d.EventAdd()
    log("OK : %d points, %d polygones." % (vase.GetPointCount(), vase.GetPolygonCount()))
    log("Export : selectionner le null MODEL puis Fichier > Exporter (FBX / OBJ / USD).")


if __name__ == '__main__':
    main()
