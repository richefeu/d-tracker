// ------------------------------------------
// D-TRACKER GUI (raylib + raygui)
// Laboratoire 3SR, Universite Grenoble Alpes
// ------------------------------------------
//
// Application graphique pour simplifier l'utilisation de tracker.
//
//   Usage:  ./tracker-gui fichier_de_commande.trk
//
// Elle s'appuie sur libtracker.a (le coeur de tracker) pour :
//   - lire le fichier de commande .trk (via init()),
//   - lire les images TIFF/RAW (via read_image()),
//   - lancer un calcul de correlation (via les procedures particle_tracking(), ...).
// L'affichage (grains, image de fond, panneau de controle) est realise avec raylib
// et raygui.
//
// Cette premiere version est un SQUELETTE : navigation entre les fichiers DIC,
// visualisation des grains colores par NCC, fond image optionnel a resolution
// degradee, et lancement d'un calcul. L'edition manuelle des grains (drag,
// correction) viendra dans un second temps.

#include "raylib.h"

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

// Lecture/ecriture de la configuration (fichier .ini), extension raylib de raysan5.
#define RINI_IMPLEMENTATION
#include "rini.h"

// Coeur de tracker (declarations seulement ; le code est dans libtracker.a)
#include "tracker.hpp"

// Outils reutilises
#include <tiffio.h> // requis par thumbnail.hpp (methodes writeTiff) : doit preceder
#include "thumbnail.hpp"
#include "toofus/ColorTable.hpp"

// Fonte d'interface embarquee (Open Sans, OFL) - generee par 'xxd -i guifont.ttf'
#include "guifont_ttf.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// ================================================================================
// Apparence (theme clair "clean" + fonte lisible)
// ================================================================================

// Palette
static const Color COL_BG = {242, 243, 245, 255};         // fond de la zone de vue
static const Color COL_PANEL = {255, 255, 255, 255};      // fond du panneau
static const Color COL_PANEL_EDGE = {221, 225, 230, 255}; // bordures / separateurs
static const Color COL_TEXT = {43, 47, 54, 255};          // texte principal
static const Color COL_MUTED = {130, 137, 145, 255};      // texte secondaire
static const Color COL_ACCENT = {45, 125, 210, 255};      // accent (titres de section)

// Couleurs des trois positions d'un grain (vue + legende du panneau)
static const Color COL_POS_REF = {120, 110, 205, 255};  // violet : position de reference
static const Color COL_POS_PREV = {235, 150, 40, 255};  // orange : position precedente
static const Color COL_POS_CUR = {40, 175, 95, 255};    // vert   : position actuelle

static Font uiFont; // fonte d'interface (chargee au demarrage)

// ---- Parametres persistants (fichier tracker-gui.ini, lus/ecrits via rini) -------
// D'autres parametres persistes sont des globales declarees plus bas (gui_image_div,
// colorMode, colNCCMin, grainAlpha, show* = affichages par defaut) ; voir load/saveConfig.
static const char *CONFIG_FILE = "tracker-gui.ini";
static int cfgWinW = 1280, cfgWinH = 800;            // dimensions de la fenetre
static int cfgFontSize = 16;                         // taille de police des widgets
static int cfgBtnH = 22, cfgSldH = 11, cfgLblH = 17; // hauteurs des widgets du panneau
static float cfgMarkerR = 4.0f;                      // rayon de base des marqueurs (ecran)
static float cfgRotR = 11.0f;                        // rayon de l'indicateur de rotation (ecran)
static double cfgFallbackRadius = 10.0;              // rayon de grain par defaut (sans radius_pix)

// Texte avec la fonte d'interface
static void uiText(const char *s, float x, float y, float size, Color c) {
  DrawTextEx(uiFont, s, {x, y}, size, 0.0f, c);
}

// Applique un theme clair coherent a tous les widgets raygui
static void applyCleanTheme() {
  GuiSetStyle(DEFAULT, TEXT_SIZE, cfgFontSize);
  GuiSetStyle(DEFAULT, TEXT_SPACING, 0);
  GuiSetStyle(DEFAULT, BORDER_WIDTH, 1);
  GuiSetStyle(DEFAULT, BACKGROUND_COLOR, 0xffffffff);
  GuiSetStyle(DEFAULT, LINE_COLOR, 0xdde1e6ff);

  GuiSetStyle(DEFAULT, BORDER_COLOR_NORMAL, 0xdde1e6ff);
  GuiSetStyle(DEFAULT, BASE_COLOR_NORMAL, 0xffffffff);
  GuiSetStyle(DEFAULT, TEXT_COLOR_NORMAL, 0x2b2f36ff);

  GuiSetStyle(DEFAULT, BORDER_COLOR_FOCUSED, 0x2d7dd2ff);
  GuiSetStyle(DEFAULT, BASE_COLOR_FOCUSED, 0xeaf2fbff);
  GuiSetStyle(DEFAULT, TEXT_COLOR_FOCUSED, 0x1b4a7aff);

  GuiSetStyle(DEFAULT, BORDER_COLOR_PRESSED, 0x2d7dd2ff);
  GuiSetStyle(DEFAULT, BASE_COLOR_PRESSED, 0x2d7dd2ff);
  GuiSetStyle(DEFAULT, TEXT_COLOR_PRESSED, 0xffffffff);

  GuiSetStyle(DEFAULT, BORDER_COLOR_DISABLED, 0xe6e9edff);
  GuiSetStyle(DEFAULT, BASE_COLOR_DISABLED, 0xf2f3f5ff);
  GuiSetStyle(DEFAULT, TEXT_COLOR_DISABLED, 0x9aa0a6ff);

  GuiSetStyle(TOGGLE, GROUP_PADDING, 2); // espacement de la grille de champs (COLORING)
  GuiSetStyle(DEFAULT, TEXT_ALIGNMENT, TEXT_ALIGN_CENTER);
  GuiSetStyle(LABEL, TEXT_ALIGNMENT, TEXT_ALIGN_LEFT);
  GuiSetStyle(SLIDER, TEXT_ALIGNMENT, TEXT_ALIGN_RIGHT);
}

// ================================================================================
// Modele d'affichage
// ================================================================================

// Un grain tel que lu dans un fichier dic_out_N.txt (13 colonnes).
struct GuiGrain {
  double refcoord_xpix, refcoord_ypix;
  double refrot, radius_pix;
  double dx, dy, drot;
  double upix, vpix, rot_inc;
  double NCC, NCC_rescue, NCC_subpix;
};

// Champs scalaires d'un grain disponibles pour la coloration. `member == nullptr`
// designe "grey" (aplat, aucun champ). `ncc == true` : valeur dans [0,1], la borne
// basse est pilotee par le slider colNCCMin ; sinon la plage est min/max automatique
// sur les grains affiches. Pour ajouter un champ, il suffit d'une ligne ici.
struct ColorField {
  const char *name;
  double GuiGrain::*member;
  bool ncc;
};
static const ColorField COLOR_FIELDS[] = {
    {"grey", nullptr, false},          {"NCC", &GuiGrain::NCC, true},
    {"NCC_rescue", &GuiGrain::NCC_rescue, true}, {"NCC_subpix", &GuiGrain::NCC_subpix, true},
    {"dx", &GuiGrain::dx, false},      {"dy", &GuiGrain::dy, false},
    {"drot", &GuiGrain::drot, false},  {"upix", &GuiGrain::upix, false},
    {"vpix", &GuiGrain::vpix, false},  {"rot_inc", &GuiGrain::rot_inc, false},
};
static const int COLOR_FIELD_COUNT = (int)(sizeof(COLOR_FIELDS) / sizeof(COLOR_FIELDS[0]));

static std::vector<GuiGrain> g_grains;
static int g_dicNum = 0;
static bool g_fromDicFile = false; // true: grains issus d'un dic_out ; false: du vecteur grain global

static int colorMode = 1;            // index dans COLOR_FIELDS (1 = NCC par defaut)
static float colNCCMin = 0.7f;       // borne basse de la table pour les champs NCC
static float grainAlpha = 1.0f;      // transparence du remplissage des disques (0=transparent, 1=opaque)
static int iSelected = -1;           // grain selectionne (clic) ; -1 = aucun
static bool showGrains = true;       // false: disques remplaces par une croix de taille fixe
static bool showTrajectory = false;  // afficher position ref + precedente (+ trace vers l'actuelle)
static bool trajSelectedOnly = false;// limiter la trace au grain selectionne
static bool showPattern = false;     // superposer le pattern correle sur chaque grain
static bool showSearchZones = false; // afficher les zones de recherche du grain selectionne

// Visualiseur (lecture seule) du fichier de commande .trk
static bool showTrk = false;
static std::string g_trkPath;
static std::string g_trkText;
static float g_trkScroll = 0.0f;

static bool showBackground = false;
static int gui_image_div = 8; // sous-echantillonnage du fond (resolution degradee)
static Texture2D bgTexture = {};
static bool bgLoaded = false;
static int bgLoadedNum = -1;
static int bgLoadedDiv = -1;

static const int PANEL_W = 270; // largeur du panneau de controle (a droite)

// Rayon d'affichage/selection par defaut, calcule quand radius_pix n'est pas donne.
static double g_fallbackRadius = 10.0;
static const float PATTERN_FALLBACK_SIDE = 15.0f; // carre 15x15 si aucun pattern defini

// Bandeau de resultat affiche apres un Run (liste des fichiers produits, statut).
static std::vector<std::string> g_resultLines;
static double g_resultUntil = 0.0; // GetTime() jusqu'auquel afficher le bandeau

// Charge les parametres depuis l'INI ; toute cle absente conserve la valeur par defaut.
static void loadConfig() {
  rini_data ini = rini_load(CONFIG_FILE); // alloue meme si le fichier n'existe pas
  cfgWinW = rini_get_value_fallback(ini, "window_width", cfgWinW);
  cfgWinH = rini_get_value_fallback(ini, "window_height", cfgWinH);
  cfgFontSize = rini_get_value_fallback(ini, "widget_font_size", cfgFontSize);
  cfgBtnH = rini_get_value_fallback(ini, "widget_button_height", cfgBtnH);
  cfgSldH = rini_get_value_fallback(ini, "widget_slider_height", cfgSldH);
  cfgLblH = rini_get_value_fallback(ini, "widget_label_height", cfgLblH);
  cfgMarkerR = (float)rini_get_value_fallback(ini, "marker_radius", (int)cfgMarkerR);
  cfgRotR = (float)rini_get_value_fallback(ini, "marker_rotation_radius", (int)cfgRotR);
  cfgFallbackRadius = rini_get_value_fallback(ini, "marker_fallback_radius", (int)cfgFallbackRadius);
  gui_image_div = rini_get_value_fallback(ini, "image_div", gui_image_div);
  colorMode = rini_get_value_fallback(ini, "color_field", colorMode);
  showGrains = rini_get_value_fallback(ini, "show_grains", showGrains) != 0;
  showTrajectory = rini_get_value_fallback(ini, "show_positions", showTrajectory) != 0;
  trajSelectedOnly = rini_get_value_fallback(ini, "positions_selected_only", trajSelectedOnly) != 0;
  showPattern = rini_get_value_fallback(ini, "show_pattern", showPattern) != 0;
  showSearchZones = rini_get_value_fallback(ini, "show_search_zones", showSearchZones) != 0;
  showBackground = rini_get_value_fallback(ini, "show_background", showBackground) != 0;
  colNCCMin = (float)atof(rini_get_value_text_fallback(ini, "ncc_min", TextFormat("%.3f", colNCCMin)));
  grainAlpha = (float)atof(rini_get_value_text_fallback(ini, "fill_opacity", TextFormat("%.3f", grainAlpha)));
  rini_unload(&ini);
}

// Enregistre les parametres courants dans l'INI (cree ou met a jour le fichier).
static void saveConfig() {
  rini_data ini = rini_load(CONFIG_FILE); // repart de l'existant : conserve l'ordre/les cles
  rini_set_value(&ini, "window_width", cfgWinW, "Largeur de la fenetre (px)");
  rini_set_value(&ini, "window_height", cfgWinH, "Hauteur de la fenetre (px)");
  rini_set_value(&ini, "widget_font_size", cfgFontSize, "Taille de police des widgets");
  rini_set_value(&ini, "widget_button_height", cfgBtnH, "Hauteur des boutons/toggles (px)");
  rini_set_value(&ini, "widget_slider_height", cfgSldH, "Hauteur des sliders (px)");
  rini_set_value(&ini, "widget_label_height", cfgLblH, "Pas vertical des lignes de texte (px)");
  rini_set_value(&ini, "marker_radius", (int)cfgMarkerR, "Rayon de base des marqueurs (ecran)");
  rini_set_value(&ini, "marker_rotation_radius", (int)cfgRotR, "Rayon de l'indicateur de rotation");
  rini_set_value(&ini, "marker_fallback_radius", (int)cfgFallbackRadius, "Rayon de grain par defaut");
  rini_set_value(&ini, "image_div", gui_image_div, "Sous-echantillonnage de l'image de fond");
  rini_set_value(&ini, "color_field", colorMode, "Champ de coloration (index COLOR_FIELDS)");
  rini_set_value(&ini, "show_grains", showGrains, "Afficher les disques (0/1)");
  rini_set_value(&ini, "show_positions", showTrajectory, "Afficher ref/prec/actuelle (0/1)");
  rini_set_value(&ini, "positions_selected_only", trajSelectedOnly, "Positions: grain selectionne seul (0/1)");
  rini_set_value(&ini, "show_pattern", showPattern, "Afficher le pattern (0/1)");
  rini_set_value(&ini, "show_search_zones", showSearchZones, "Afficher les zones de recherche (0/1)");
  rini_set_value(&ini, "show_background", showBackground, "Afficher l'image de fond (0/1)");
  rini_set_value_text(&ini, "ncc_min", TextFormat("%.3f", colNCCMin), "Borne basse NCC [0..1]");
  rini_set_value_text(&ini, "fill_opacity", TextFormat("%.3f", grainAlpha), "Opacite du remplissage [0..1]");
  rini_save(ini, CONFIG_FILE);
  rini_unload(&ini);
  TraceLog(LOG_INFO, "Settings saved to %s", CONFIG_FILE);
}

// ================================================================================
// Rayon effectif d'un point (affichage et selection)
// ================================================================================

// Calcule un rayon par defaut a partir de la taille de l'image (si connue) ou de la
// bounding box des points, et du nombre de points suivis.
static void computeFallbackRadius() {
  double area = 0.0;
  if (dimx > 0 && dimy > 0) {
    area = (double)dimx * (double)dimy;
  } else if (!g_grains.empty()) {
    double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;
    for (const auto &g : g_grains) {
      double x = g.refcoord_xpix + g.dx, y = g.refcoord_ypix + g.dy;
      if (x < xmin)
        xmin = x;
      if (x > xmax)
        xmax = x;
      if (y < ymin)
        ymin = y;
      if (y > ymax)
        ymax = y;
    }
    if (xmax > xmin && ymax > ymin)
      area = (xmax - xmin) * (ymax - ymin);
  }
  int n = (int)g_grains.size();
  if (area > 0.0 && n > 0)
    g_fallbackRadius = 0.35 * std::sqrt(area / (double)n); // ~ fraction de l'espacement moyen
  else
    g_fallbackRadius = cfgFallbackRadius;
}

// Rayon a utiliser pour le point i :
//   - radius_pix s'il est fourni (> 0) ;
//   - sinon, meme surface que le pattern s'il existe : r = sqrt(N_pix / pi) ;
//   - sinon, rayon par defaut calcule (image / nombre de points).
static double effectiveRadius(int i) {
  double r = g_grains[i].radius_pix;
  if (r > 0.0)
    return r;
  if (i < (int)grain.size() && !grain[i].pattern.empty())
    return std::sqrt((double)grain[i].pattern.size() / M_PI);
  return g_fallbackRadius;
}

// ================================================================================
// Lecture d'un fichier DIC (porte depuis 1g2e_tools/seeDIC/dicProcess.cpp)
// ================================================================================

static bool fileExists(const char *name) {
  std::ifstream f(name);
  return f.good();
}

// Construit le nom du fichier DIC a partir du format dic_name (printf-style)
static void dicFileName(int num, char *out, size_t n) {
  // dic_name est par defaut quelque chose comme "dic_out_%d.txt"
  if (dic_name[0] != '\0')
    snprintf(out, n, dic_name, num);
  else
    snprintf(out, n, "dic_out_%d.txt", num);
}

static bool readDIC(int num) {
  char fname[512];
  dicFileName(num, fname, sizeof(fname));
  if (!fileExists(fname)) {
    TraceLog(LOG_WARNING, "DIC file does not exist: %s", fname);
    return false;
  }
  std::ifstream file(fname);
  if (!file) {
    TraceLog(LOG_WARNING, "Cannot open DIC file: %s", fname);
    return false;
  }
  int ng = 0;
  file >> ng;
  g_grains.clear();
  g_grains.reserve(ng > 0 ? ng : 0);
  GuiGrain G;
  for (int i = 0; i < ng; ++i) {
    file >> G.refcoord_xpix >> G.refcoord_ypix >> G.refrot >> G.radius_pix >> G.dx >> G.dy >> G.drot >> G.upix >>
        G.vpix >> G.rot_inc >> G.NCC >> G.NCC_rescue >> G.NCC_subpix;
    g_grains.push_back(G);
  }
  g_dicNum = num;
  g_fromDicFile = true;
  // On conserve la selection (meme grain physique d'un fichier a l'autre), sauf si
  // l'indice n'existe plus dans le nouveau jeu de grains.
  if (iSelected >= (int)g_grains.size())
    iSelected = -1;
  computeFallbackRadius();
  TraceLog(LOG_INFO, "Loaded %s (%zu grains)", fname, g_grains.size());
  return true;
}

// Lecture d'un fichier de points au format "grains" (1re ligne: nombre, puis par
// ligne: x y refrot radius). C'est le format de grains_to_follow et de recentered.data.
static bool readGrainsInputFile(const char *path) {
  std::ifstream f(path);
  if (!f)
    return false;
  int n = 0;
  f >> n;
  if (!f || n <= 0)
    return false;
  g_grains.clear();
  g_grains.reserve((size_t)n);
  for (int i = 0; i < n; ++i) {
    double x, y, r, rad;
    if (!(f >> x >> y >> r >> rad))
      break;
    GuiGrain G{};
    G.refcoord_xpix = x;
    G.refcoord_ypix = y;
    G.refrot = r;
    G.radius_pix = rad;
    g_grains.push_back(G);
  }
  g_fromDicFile = false;
  if (iSelected >= (int)g_grains.size())
    iSelected = -1;
  computeFallbackRadius();
  TraceLog(LOG_INFO, "Loaded %zu grains from %s", g_grains.size(), path);
  return true;
}

// Repli quand aucun fichier dic_out n'est disponible (procedures sans suivi de type
// find_subpixel_centers, pattern_quality, ...) : on affiche directement les grains
// charges par init() dans le vecteur global `grain` (issus de grains_to_follow, etc.).
static void loadGrainsFromGlobal() {
  g_grains.clear();
  g_grains.reserve((size_t)num_grains);
  for (int i = 0; i < num_grains; ++i) {
    GuiGrain G;
    G.refcoord_xpix = grain[i].refcoord_xpix;
    G.refcoord_ypix = grain[i].refcoord_ypix;
    G.refrot = grain[i].refrot;
    G.radius_pix = grain[i].radius_pix;
    G.dx = grain[i].dx;
    G.dy = grain[i].dy;
    G.drot = grain[i].drot;
    G.upix = grain[i].upix;
    G.vpix = grain[i].vpix;
    G.rot_inc = grain[i].rot_inc;
    G.NCC = grain[i].NCC;
    G.NCC_rescue = grain[i].NCC_rescue;
    G.NCC_subpix = grain[i].NCC_subpix;
    g_grains.push_back(G);
  }
  g_fromDicFile = false;
  if (iSelected >= (int)g_grains.size())
    iSelected = -1; // conserve la selection si possible
  computeFallbackRadius();
  TraceLog(LOG_INFO, "Loaded %d grains from the command file (no dic_out)", num_grains);
}

// ================================================================================
// Fond image (chargee via libtracker.a, sous-echantillonnee via thumbnail)
// ================================================================================

static void unloadBackground() {
  if (bgLoaded) {
    UnloadTexture(bgTexture);
    bgLoaded = false;
    bgLoadedNum = -1;
    bgLoadedDiv = -1;
  }
}

// Charge l'image numero `num` et en construit une texture sous-echantillonnee.
static void loadBackground(int num) {
  // Verifie d'abord que le fichier image existe (read_image fait exit(1) sinon).
  char name[512];
  if (image_name[0] == '\0') {
    TraceLog(LOG_WARNING, "No image_name defined in the command file");
    return;
  }
  snprintf(name, sizeof(name), image_name, num);
  if (!fileExists(name)) {
    TraceLog(LOG_WARNING, "Background image does not exist: %s", name);
    return;
  }

  // Lecture de l'image dans le tableau global image[0] (via libtracker.a)
  read_image(0, num, /*first_time=*/true);
  if (dimx <= 0 || dimy <= 0) {
    TraceLog(LOG_WARNING, "Invalid image size for %s", name);
    return;
  }

  if (gui_image_div < 1)
    gui_image_div = 1;
  thumbnail th(image, dimx, dimy, 0, gui_image_div);

  // Conversion thumbnail (uint16, colonne-major) -> Image raylib RGBA (row-major)
  int W = th.imW;
  int H = th.imH;
  unsigned char *pixels = (unsigned char *)MemAlloc((unsigned int)(W * H * 4));
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      unsigned char v = (unsigned char)((th.imThumb[x][y] * 255) / 65535);
      int idx = (y * W + x) * 4;
      pixels[idx + 0] = v;
      pixels[idx + 1] = v;
      pixels[idx + 2] = v;
      pixels[idx + 3] = 255;
    }
  }
  Image img;
  img.data = pixels;
  img.width = W;
  img.height = H;
  img.mipmaps = 1;
  img.format = PIXELFORMAT_UNCOMPRESSED_R8G8B8A8;

  unloadBackground();
  bgTexture = LoadTextureFromImage(img);
  UnloadImage(img); // libere pixels
  bgLoaded = true;
  bgLoadedNum = num;
  bgLoadedDiv = gui_image_div;
  TraceLog(LOG_INFO, "Background loaded (%dx%d, div=%d)", W, H, gui_image_div);
}

// ================================================================================
// Vue / camera
// ================================================================================

// Recadre la vue sur l'ensemble des grains.
static void fitView(Camera2D &cam, int screenW, int screenH) {
  if (g_grains.empty())
    return;
  double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;
  for (size_t i = 0; i < g_grains.size(); ++i) {
    const GuiGrain &g = g_grains[i];
    double xc = g.refcoord_xpix + g.dx;
    double yc = g.refcoord_ypix + g.dy;
    double R = effectiveRadius((int)i);
    if (xc - R < xmin)
      xmin = xc - R;
    if (xc + R > xmax)
      xmax = xc + R;
    if (yc - R < ymin)
      ymin = yc - R;
    if (yc + R > ymax)
      ymax = yc + R;
  }
  double w = xmax - xmin;
  double h = ymax - ymin;
  if (w <= 0)
    w = 1;
  if (h <= 0)
    h = 1;
  int viewW = screenW - PANEL_W;
  double zx = viewW / (w * 1.1);
  double zy = screenH / (h * 1.1);
  cam.zoom = (float)((zx < zy) ? zx : zy);
  cam.target = {(float)(0.5 * (xmin + xmax)), (float)(0.5 * (ymin + ymax))};
  cam.offset = {(float)(viewW * 0.5f), (float)(screenH * 0.5f)};
}

// ================================================================================
// Dessin
// ================================================================================

static Color ncc2color(ColorTable &CT, float value) {
  colorRGBA c;
  CT.getRGB(value, &c);
  return Color{(unsigned char)c.r, (unsigned char)c.g, (unsigned char)c.b, 255};
}

// Plage [lo, hi] de la table de couleurs pour le champ courant : [colNCCMin, 1] pour
// un champ NCC, sinon min/max automatique sur les grains affiches.
static void colorRange(const ColorField &cf, float &lo, float &hi) {
  if (cf.ncc || !cf.member) {
    lo = colNCCMin;
    hi = 1.0f;
    return;
  }
  double mn = 1e30, mx = -1e30;
  for (const auto &g : g_grains) {
    double v = g.*(cf.member);
    if (v < mn)
      mn = v;
    if (v > mx)
      mx = v;
  }
  if (mx <= mn)
    mx = mn + 1.0; // plage non nulle pour eviter une division par zero
  lo = (float)mn;
  hi = (float)mx;
}

static void drawGrains(float zoom) {
  const ColorField &cf = COLOR_FIELDS[colorMode];
  float lo, hi;
  colorRange(cf, lo, hi);

  ColorTable CT;
  CT.setSize(16);
  CT.setSwap(cf.ncc); // pour NCC : couleur "chaude" cote 1 (bonne correlation)
  CT.setMinMax(lo, hi);
  CT.Rebuild();

  // Demi-taille et epaisseur de la croix (mode "grains masques"), constantes a l'ecran.
  float cross = (zoom > 0.0f) ? cfgMarkerR / zoom : cfgMarkerR;
  float crossTh = (zoom > 0.0f) ? 1.5f / zoom : 1.5f;

  for (size_t i = 0; i < g_grains.size(); ++i) {
    const GuiGrain &g = g_grains[i];
    float xc = (float)(g.refcoord_xpix + g.dx);
    float yc = (float)(g.refcoord_ypix + g.dy);
    // radius_pix s'il est donne, sinon rayon effectif (surface du pattern ou calcule)
    float R = (float)effectiveRadius((int)i);
    Vector2 c = {xc, yc};

    if (R == 1.0f) { // Corners -> croix rouge
      DrawLineEx({xc - 40, yc - 40}, {xc + 40, yc + 40}, 2.0f, RED);
      DrawLineEx({xc - 40, yc + 40}, {xc + 40, yc - 40}, 2.0f, RED);
      continue;
    } else if (R == 2.1f) { // scaling points -> croix bleue
      DrawLineEx({xc, yc - 40}, {xc, yc + 40}, 2.0f, BLUE);
      DrawLineEx({xc - 40, yc}, {xc + 40, yc}, 2.0f, BLUE);
      continue;
    } else if (R == 2.2f) { // fixed points -> point vert
      DrawCircleV(c, 6.0f, GREEN);
      continue;
    }

    Color fill;
    if (cf.member)
      fill = ncc2color(CT, (float)(g.*(cf.member)));
    else
      fill = Color{230, 230, 230, 255}; // "grey" : aplat neutre

    if (!showGrains) {
      // Grains masques : croix de taille fixe (ecran) coloree par le champ courant.
      fill.a = 255;
      DrawLineEx({xc - cross, yc - cross}, {xc + cross, yc + cross}, crossTh, fill);
      DrawLineEx({xc - cross, yc + cross}, {xc + cross, yc - cross}, crossTh, fill);
      if ((int)i == iSelected)
        DrawCircleLinesV(c, cross + 3.0f * crossTh, ORANGE);
      continue;
    }

    fill.a = (unsigned char)(grainAlpha * 255.0f); // remplissage transparent reglable
    DrawCircleV(c, R, fill);
    DrawCircleLinesV(c, R, BLACK); // contour toujours opaque

    if ((int)i == iSelected) {
      DrawRing(c, R, R + 3.0f, 0, 360, 32, ORANGE);
    }
  }
}

// Colorbar verticale en haut a droite de la vue, sauf si le champ "grey" est choisi.
// Affiche le degrade du champ courant et ses bornes (haut = max, bas = min).
static void drawColorBar(int viewW) {
  const ColorField &cf = COLOR_FIELDS[colorMode];
  if (!cf.member)
    return; // aucune colorbar pour "grey"

  float lo, hi;
  colorRange(cf, lo, hi);
  ColorTable CT;
  CT.setSize(16);
  CT.setSwap(cf.ncc);
  CT.setMinMax(lo, hi);
  CT.Rebuild();

  const float bw = 18.0f, bh = 170.0f;
  const float bx = viewW - 24.0f - bw; // ancrage a droite de la vue
  const float by = 38.0f;

  // Carte de fond pour lisibilite
  DrawRectangleRounded(Rectangle{bx - 10, by - 26, bw + 20, bh + 36}, 0.12f, 8, Fade(COL_PANEL, 0.94f));

  // Degrade : bandes du haut (max) vers le bas (min)
  const int N = 64;
  for (int k = 0; k < N; ++k) {
    float t = (k + 0.5f) / N;             // 0 en haut, 1 en bas
    float v = lo + (hi - lo) * (1.0f - t); // haut = hi
    Color col = ncc2color(CT, v);
    DrawRectangle((int)bx, (int)(by + bh * k / N), (int)bw, (int)(bh / N + 1.5f), col);
  }
  DrawRectangleLinesEx(Rectangle{bx, by, bw, bh}, 1.0f, COL_PANEL_EDGE);

  // Titre (nom du champ) au-dessus, aligne a droite du bandeau
  Vector2 ns = MeasureTextEx(uiFont, cf.name, 15, 0);
  uiText(cf.name, bx + bw - ns.x, by - 22, 15, COL_TEXT);

  // Bornes max (haut) et min (bas), a gauche du bandeau, alignees a droite
  auto rightLabel = [&](const char *s, float yy) {
    Vector2 ts = MeasureTextEx(uiFont, s, 14, 0);
    uiText(s, bx - 6 - ts.x, yy, 14, COL_MUTED);
  };
  rightLabel(TextFormat("%.3g", hi), by - 4);
  rightLabel(TextFormat("%.3g", lo), by + bh - 12);
}

// Positions ref / precedente / actuelle d'un grain :
//   ref       = (refcoord_xpix, refcoord_ypix)
//   actuelle  = ref + (dx, dy)
//   precedente= actuelle - (upix, vpix)  [upix/vpix = increment du dernier pas DIC]
static void grainPositions(const GuiGrain &g, Vector2 &ref, Vector2 &prev, Vector2 &cur) {
  ref = {(float)g.refcoord_xpix, (float)g.refcoord_ypix};
  cur = {(float)(g.refcoord_xpix + g.dx), (float)(g.refcoord_ypix + g.dy)};
  prev = {(float)(cur.x - g.upix), (float)(cur.y - g.vpix)};
}

// Trace, pour tous les grains (ou seulement le selectionne), la position de reference
// (violet) et la position precedente (orange) reliees a la position actuelle (vert).
// Les marqueurs ont une taille constante a l'ecran. A dessiner dans BeginMode2D.
static void drawTrajectories(float zoom) {
  if (!showTrajectory)
    return;
  float r = (zoom > 0.0f) ? cfgMarkerR / zoom : cfgMarkerR; // rayon des marqueurs (ecran)
  float th = (zoom > 0.0f) ? 1.5f / zoom : 1.5f;            // epaisseur des traits

  float rotR = (zoom > 0.0f) ? cfgRotR / zoom : cfgRotR; // rayon de l'indicateur de rotation

  for (size_t i = 0; i < g_grains.size(); ++i) {
    if (trajSelectedOnly && (int)i != iSelected)
      continue;
    const GuiGrain &g = g_grains[i];
    Vector2 ref, prev, cur;
    grainPositions(g, ref, prev, cur);

    // Trace ref -> precedente -> actuelle
    DrawLineEx(ref, prev, th, Fade(COL_POS_REF, 0.6f));
    DrawLineEx(prev, cur, th, COL_POS_PREV);

    // Marqueurs ref (cercle violet) et precedente (cercle orange)
    DrawCircleLinesV(ref, r, COL_POS_REF);
    DrawCircleV(ref, r * 0.35f, COL_POS_REF);
    DrawCircleLinesV(prev, r, COL_POS_PREV);

    // ---- Position actuelle : indicateur de rotation ----
    // Grand cercle + croix en X a l'orientation de reference (refrot, croix violette non
    // tournee) et a l'orientation actuelle (refrot + drot, croix verte tournee) ; le
    // secteur rempli entre les deux visualise l'angle de rotation (drot).
    auto drawX = [&](float ang, float thick, Color col) {
      for (int k = 0; k < 2; ++k) {
        float a = ang + (float)(M_PI / 4.0) + k * (float)(M_PI / 2.0); // bras a 45 et 135
        Vector2 d = {cosf(a) * rotR, sinf(a) * rotR};
        DrawLineEx({cur.x - d.x, cur.y - d.y}, {cur.x + d.x, cur.y + d.y}, thick, col);
      }
    };
    float a0 = (float)((g.refrot + M_PI / 4.0) * RAD2DEG);
    float a1 = (float)((g.refrot + g.drot + M_PI / 4.0) * RAD2DEG);
    DrawCircleSector(cur, rotR * 0.72f, fminf(a0, a1), fmaxf(a0, a1), 24, Fade(COL_POS_CUR, 0.30f));
    DrawCircleLinesV(cur, rotR, COL_POS_CUR);
    drawX(g.refrot, th, Fade(COL_POS_REF, 0.55f));    // croix non tournee (orientation ref)
    drawX(g.refrot + g.drot, th * 1.3f, COL_POS_CUR); // croix tournee (orientation actuelle)
    DrawCircleV(cur, th * 1.3f, COL_POS_CUR);         // centre precis
  }
}

// Dessine le pattern (zone de pixels correlee) place sur chaque grain.
//
// Le pattern n'est pas stocke dans les fichiers dic_out : il provient de la commande
// 'pattern ...' du fichier .trk. Comme init() (de libtracker.a) a deja construit les
// patterns dans le vecteur global `grain`, on lit directement grain[i].pattern (liste
// d'offsets entiers relatifs au centre). L'index i correspond a g_grains[i] (meme ordre
// que save_grains). Chaque pixel du pattern est dessine comme un petit carre 1x1, place
// a la position courante du grain et tourne de (refrot + drot).
static void drawPatterns() {
  static const std::vector<relative_coord_type> emptyPat;
  for (size_t i = 0; i < g_grains.size(); ++i) {
    // Pattern reel si disponible dans le vecteur global, sinon vide -> fallback carre.
    const std::vector<relative_coord_type> &pat = (i < grain.size()) ? grain[i].pattern : emptyPat;

    const GuiGrain &g = g_grains[i];
    double cx = g.refcoord_xpix + g.dx;
    double cy = g.refcoord_ypix + g.dy;
    double ang = g.refrot + g.drot;
    double ca = cos(ang), sa = sin(ang);

    Color col = ((int)i == iSelected) ? Color{255, 140, 0, 160} : Color{0, 150, 160, 110};

    if (!pat.empty()) {
      // Vrai pattern : chaque pixel dessine comme un petit carre 1x1.
      for (const relative_coord_type &p : pat) {
        double rx = p.dx * ca - p.dy * sa;
        double ry = p.dx * sa + p.dy * ca;
        DrawRectangleV({(float)(cx + rx - 0.5), (float)(cy + ry - 0.5)}, {1.0f, 1.0f}, col);
      }
    } else {
      // Aucun pattern defini dans le .trk : carre indicatif 15x15 (tourne) en contour.
      float h = PATTERN_FALLBACK_SIDE * 0.5f;
      Vector2 corners[4];
      float lx[4] = {-h, h, h, -h};
      float ly[4] = {-h, -h, h, h};
      for (int k = 0; k < 4; ++k) {
        corners[k].x = (float)(cx + lx[k] * ca - ly[k] * sa);
        corners[k].y = (float)(cy + lx[k] * sa + ly[k] * ca);
      }
      for (int k = 0; k < 4; ++k)
        DrawLineV(corners[k], corners[(k + 1) % 4], col);
    }
  }
}

// Dessine, autour du grain selectionne, les rectangles des zones de recherche
// (normale / rescue / super-rescue). Les dimensions proviennent des globales
// search_zone* (chargees par init() depuis le .trk). A dessiner dans BeginMode2D.
// `zoom` sert a garder une epaisseur de trait ~constante a l'ecran.
static void drawSearchZones(float zoom) {
  if (iSelected < 0 || iSelected >= (int)g_grains.size())
    return;
  const GuiGrain &g = g_grains[iSelected];
  float cx = (float)(g.refcoord_xpix + g.dx);
  float cy = (float)(g.refcoord_ypix + g.dy);
  float t = (zoom > 0.0f) ? 1.5f / zoom : 1.5f;

  // left/right vers x-/x+, up/down vers y-/y+ (cf. struct search_zone_type)
  auto box = [&](const search_zone_type &z) {
    return Rectangle{cx - (float)z.left, cy - (float)z.up, (float)(z.left + z.right), (float)(z.up + z.down)};
  };
  // De la plus grande a la plus petite pour que les petites restent visibles
  DrawRectangleLinesEx(box(search_zone_super_rescue), t, Color{220, 60, 60, 255}); // rouge
  DrawRectangleLinesEx(box(search_zone_rescue), t, Color{235, 140, 30, 255});      // orange
  DrawRectangleLinesEx(box(search_zone), t, Color{40, 170, 80, 255});              // vert
}

// Selection du grain sous le curseur (coordonnees monde).
static void pickGrain(Vector2 world) {
  iSelected = -1;
  for (size_t i = 0; i < g_grains.size(); ++i) {
    const GuiGrain &g = g_grains[i];
    double xc = g.refcoord_xpix + g.dx;
    double yc = g.refcoord_ypix + g.dy;
    double R = effectiveRadius((int)i); // meme rayon que l'affichage (radius_pix ou calcule)
    double ddx = xc - world.x;
    double ddy = yc - world.y;
    if (ddx * ddx + ddy * ddy < R * R) {
      iSelected = (int)i;
      return;
    }
  }
}

// Panneau flottant listant TOUTES les informations disponibles sur le grain
// selectionne. Dessine en coordonnees ecran (apres EndMode2D).
static void drawInspector() {
  if (iSelected < 0 || iSelected >= (int)g_grains.size())
    return;
  const GuiGrain &g = g_grains[iSelected];
  double cx = g.refcoord_xpix + g.dx;
  double cy = g.refcoord_ypix + g.dy;
  bool hasPattern = (iSelected < (int)grain.size());

  const float ix = 16.0f, iw = 310.0f, lineH = 19.0f, top = 58.0f;
  int nL = 12 + (g_fromDicFile ? 1 : 0) + (hasPattern ? 1 : 0);
  float cardH = 30.0f + nL * lineH + 8.0f;
  DrawRectangleRounded(Rectangle{ix - 6, top - 6, iw, cardH}, 0.05f, 8, Fade(COL_PANEL, 0.95f));

  float y = top;
  uiText("SELECTED GRAIN", ix, y, 14, COL_ACCENT);
  y += 24;

  // label (muted) a gauche, valeur (texte) en colonne fixe
  auto row = [&](const char *label, const char *value) {
    uiText(label, ix, y, 15, COL_MUTED);
    uiText(value, ix + 118, y, 15, COL_TEXT);
    y += lineH;
  };

  row("index", TextFormat("%d", iSelected));
  if (g_fromDicFile)
    row("file line", TextFormat("%d  (dic_out_%d)", 2 + iSelected, g_dicNum));
  row("refcoord", TextFormat("%.0f , %.0f px", g.refcoord_xpix, g.refcoord_ypix));
  row("radius", TextFormat("%.4f px", g.radius_pix));
  row("refrot", TextFormat("%.6f rad", g.refrot));
  row("position", TextFormat("%.3f , %.3f px", cx, cy));
  row("dx , dy", TextFormat("%.4f , %.4f", g.dx, g.dy));
  row("drot", TextFormat("%.6f rad", g.drot));
  row("upix , vpix", TextFormat("%.4f , %.4f", g.upix, g.vpix));
  row("rot_inc", TextFormat("%.6f rad", g.rot_inc));
  row("NCC", TextFormat("%.5f", g.NCC));
  row("NCC_rescue", TextFormat("%.5f", g.NCC_rescue));
  row("NCC_subpix", TextFormat("%.5f", g.NCC_subpix));
  if (hasPattern)
    row("pattern pts", TextFormat("%zu", grain[iSelected].pattern.size()));
}

// Legende (espace ecran) des zones de recherche, avec leurs dimensions et parametres
// de rotation. Affichee en bas a gauche quand un grain est selectionne.
static void drawSearchZonesLegend(int sh) {
  if (iSelected < 0 || iSelected >= (int)g_grains.size())
    return;

  const float lx = 16.0f, lw = 360.0f, lineH = 19.0f;
  const int nL = 4;
  float cardH = 30.0f + nL * lineH + 8.0f;
  float top = sh - 40.0f - cardH;
  DrawRectangleRounded(Rectangle{lx - 6, top - 6, lw, cardH}, 0.05f, 8, Fade(COL_PANEL, 0.95f));

  float y = top;
  uiText("SEARCH ZONES (selected grain)", lx, y, 14, COL_ACCENT);
  y += 24;

  auto zrow = [&](Color sw, const char *name, const search_zone_type &z) {
    DrawRectangle((int)lx, (int)(y + 2), 12, 12, sw);
    uiText(name, lx + 20, y, 15, COL_TEXT);
    uiText(TextFormat("L/R/U/D %d/%d/%d/%d   rot %dx%.4f", z.left, z.right, z.up, z.down, z.num_rot, z.inc_rot),
           lx + 110, y, 15, COL_MUTED);
    y += lineH;
  };
  zrow(Color{40, 170, 80, 255}, "normal", search_zone);
  zrow(Color{235, 140, 30, 255}, "rescue", search_zone_rescue);
  zrow(Color{220, 60, 60, 255}, "super", search_zone_super_rescue);
}

// Visualiseur (lecture seule) du fichier de commande .trk : overlay scrollable.
static void drawTrkViewer(int viewW, int sh) {
  if (!showTrk)
    return;

  Rectangle area = {40, 40, (float)viewW - 80, (float)sh - 80};
  DrawRectangleRounded(area, 0.02f, 10, Fade(COL_PANEL, 0.98f));
  DrawRectangleRoundedLinesEx(area, 0.02f, 10, 1.0f, COL_PANEL_EDGE);

  uiText(TextFormat("Command file:  %s", g_trkPath.c_str()), area.x + 16, area.y + 12, 16, COL_ACCENT);
  uiText("read-only  -  wheel to scroll  -  T or Esc to close", area.x + 16, area.y + 34, 13, COL_MUTED);

  Rectangle textArea = {area.x + 16, area.y + 58, area.width - 32, area.height - 74};
  float contentH = MeasureTextEx(uiFont, g_trkText.c_str(), 15, 0).y;
  float maxScroll = (contentH > textArea.height) ? (contentH - textArea.height) : 0.0f;
  if (g_trkScroll < 0.0f)
    g_trkScroll = 0.0f;
  if (g_trkScroll > maxScroll)
    g_trkScroll = maxScroll;

  BeginScissorMode((int)textArea.x, (int)textArea.y, (int)textArea.width, (int)textArea.height);
  DrawTextEx(uiFont, g_trkText.c_str(), {textArea.x, textArea.y - g_trkScroll}, 15, 0, COL_TEXT);
  EndScissorMode();
}

// ================================================================================
// Lancement d'un calcul (memes procedures que run.cpp, via libtracker.a)
// ================================================================================

// Libelle du bouton Run adapte a la procedure courante.
static const char *procedureLabel() {
  if (procedure == "particle_tracking")
    return "Run tracking";
  if (procedure == "particle_tracking_assisted_corrections")
    return "Run corrections";
  if (procedure == "correction_distortion")
    return "Correct distortion";
  if (procedure == "correction_distortion_grid")
    return "Correct distortion";
  if (procedure == "find_subpixel_centers")
    return "Find centers";
  if (procedure == "pattern_quality")
    return "Pattern quality";
  if (procedure == "post_process")
    return "Post-process";
  if (procedure == "visu_process")
    return "Generate visu";
  if (procedure == "gray_level_analysis")
    return "Gray analysis";
  return "Run";
}

// Instantane des fichiers du repertoire courant (nom -> date de modification),
// pour detecter ce qu'une procedure a produit ou modifie.
static std::map<std::string, std::filesystem::file_time_type> snapshotDir() {
  std::map<std::string, std::filesystem::file_time_type> m;
  std::error_code ec;
  for (auto it = std::filesystem::directory_iterator(".", ec); !ec && it != std::filesystem::directory_iterator();
       it.increment(ec)) {
    if (it->is_regular_file(ec))
      m[it->path().filename().string()] = it->last_write_time(ec);
  }
  return m;
}

static void dispatchProcedure() {
  if (procedure == "particle_tracking_assisted_corrections") {
    particle_tracking_assisted_corrections();
  } else if (procedure == "correction_distortion") {
    correction_distortion();
  } else if (procedure == "correction_distortion_grid") {
    correction_distortion_grid();
  } else if (procedure == "find_subpixel_centers") {
    find_subpixel_centers();
  } else if (procedure == "pattern_quality") {
    pattern_quality();
  } else if (procedure == "post_process") {
    post_process();
  } else if (procedure == "visu_process") {
    visu_process();
  } else if (procedure == "gray_level_analysis") {
    gray_level_analysis();
  } else {
    particle_tracking();
  }
}

static void runCorrelation() {
  TraceLog(LOG_INFO, "=== Running procedure: %s ===", procedure.c_str());
  auto before = snapshotDir();

  dispatchProcedure();

  TraceLog(LOG_INFO, "=== Procedure finished ===");

  // Fichiers produits ou modifies pendant l'execution.
  auto after = snapshotDir();
  std::vector<std::string> produced;
  for (const auto &kv : after) {
    auto it = before.find(kv.first);
    if (it == before.end() || kv.second > it->second)
      produced.push_back(kv.first);
  }
  std::sort(produced.begin(), produced.end());

  // Rechargement automatique de ce qui est spatial :
  //   1) dic_out du numero courant si present ;
  //   2) sinon recentered.data (find_subpixel_centers) ;
  //   3) sinon les grains du vecteur global (mis a jour par la procedure).
  if (!readDIC(g_dicNum)) {
    if (fileExists("recentered.data"))
      readGrainsInputFile("recentered.data");
    else
      loadGrainsFromGlobal();
  }
  if (showBackground)
    loadBackground(g_dicNum);

  // Bandeau recapitulatif (affiche quelques secondes).
  g_resultLines.clear();
  g_resultLines.push_back(TextFormat("%s  -  finished", procedure.c_str()));
  if (produced.empty()) {
    g_resultLines.push_back("no file produced or modified");
  } else {
    g_resultLines.push_back(TextFormat("produced / updated %d file(s):", (int)produced.size()));
    std::string line;
    int shown = 0;
    for (const auto &f : produced) {
      if (shown >= 8) {
        line += TextFormat("  (+%d more)", (int)produced.size() - shown);
        break;
      }
      if (!line.empty())
        line += "   ";
      line += f;
      ++shown;
    }
    g_resultLines.push_back(line);
  }
  g_resultUntil = GetTime() + 12.0;
}

// Bandeau de resultat (espace ecran), centre en bas de la zone de vue.
static void drawResultBanner(int viewW, int sh) {
  if (g_resultLines.empty() || GetTime() > g_resultUntil)
    return;

  const float lineH = 20.0f;
  float w = 0.0f;
  for (const auto &s : g_resultLines) {
    float lw = MeasureTextEx(uiFont, s.c_str(), 16, 0).x;
    if (lw > w)
      w = lw;
  }
  w += 32.0f;
  if (w > viewW - 40)
    w = (float)(viewW - 40);
  float h = 20.0f + lineH * g_resultLines.size();
  float x = (viewW - w) * 0.5f;
  float y = sh - 30.0f - h;

  DrawRectangleRounded(Rectangle{x, y, w, h}, 0.12f, 8, Fade(COL_TEXT, 0.92f));
  float ty = y + 10.0f;
  bool first = true;
  for (const auto &s : g_resultLines) {
    uiText(s.c_str(), x + 16, ty, first ? 17 : 15, first ? Color{255, 255, 255, 255} : Color{210, 215, 222, 255});
    ty += lineH;
    first = false;
  }
}

// ================================================================================
// Navigation
// ================================================================================

static void gotoDIC(int num) {
  if (iinc != 0) {
    if (num < ibeg)
      num = ibeg;
    if (num > iend)
      num = iend;
  }
  if (readDIC(num)) {
    if (showBackground)
      loadBackground(num);
  }
}

// ================================================================================
// Main
// ================================================================================

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s fichier_de_commande.trk\n", argv[0]);
    return 1;
  }

  // init() (de libtracker.a) affiche le header, lit le fichier .trk et initialise
  // toutes les variables globales (image_name, dic_name, iref/ibeg/iend/iinc, ...).
  if (init(argc, argv) == 0) {
    return 0;
  }

  // Memorise le contenu brut du .trk pour le visualiseur (lecture seule).
  g_trkPath = argv[1];
  {
    std::ifstream tf(argv[1]);
    if (tf) {
      std::stringstream ss;
      ss << tf.rdbuf();
      g_trkText = ss.str();
    } else {
      g_trkText = "(impossible de relire le fichier de commande)";
    }
  }

  // Premier fichier DIC a afficher : ibeg si defini, sinon iref.
  int startNum = (iinc != 0) ? ibeg : iref;
  if (!readDIC(startNum) && !readDIC(iref) && !readDIC(1)) {
    // Aucun fichier dic_out : on affiche les grains charges par init()
    // (procedures sans suivi : find_subpixel_centers, pattern_quality, ...).
    loadGrainsFromGlobal();
    g_dicNum = (iref > 0) ? iref : 1; // numero d'image de reference (pour le fond)
  }

  loadConfig(); // parametres persistants (taille fenetre, cosmetique, affichages) avant InitWindow
  SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT);
  InitWindow(cfgWinW, cfgWinH, "D-TRACKER GUI");
  SetWindowMinSize(900, 600);
  SetExitKey(KEY_NULL); // Esc gere a la main : ferme l'overlay .trk s'il est ouvert, sinon quitte
  SetTargetFPS(60);

  // Fonte d'interface : chargee depuis la memoire (Open Sans embarquee), a une
  // taille de base elevee + mipmaps/bilineaire pour un rendu net a toute echelle.
  uiFont = LoadFontFromMemory(".ttf", guifont_ttf, (int)guifont_ttf_len, 48, NULL, 0);
  if (uiFont.texture.id != 0) {
    GenTextureMipmaps(&uiFont.texture);
    SetTextureFilter(uiFont.texture, TEXTURE_FILTER_BILINEAR);
    GuiSetFont(uiFont);
  } else {
    uiFont = GetFontDefault(); // repli improbable
  }
  applyCleanTheme();

  // Honore l'etat initial du fond (necessite le contexte GL, donc apres InitWindow).
  if (showBackground)
    loadBackground(g_dicNum);

  Camera2D cam = {{0, 0}, {0, 0}, 0.0f, 1.0f};
  fitView(cam, GetScreenWidth(), GetScreenHeight());

  while (!WindowShouldClose()) {
    int sw = GetScreenWidth();
    int sh = GetScreenHeight();
    int viewW = sw - PANEL_W;

    // ---- Input ----
    Vector2 mouse = GetMousePosition();
    bool overView = mouse.x < viewW;

    // Caracteres tapes cette frame, independants de la disposition clavier
    // (AZERTY/QWERTY) : on lit les caracteres plutot que les scancodes physiques.
    bool typed[128] = {false};
    for (int c = GetCharPressed(); c > 0; c = GetCharPressed()) {
      if (c >= 'A' && c <= 'Z')
        c += 32; // normalise en minuscule
      if (c > 0 && c < 128)
        typed[c] = true;
    }

    if (showTrk) {
      // Le visualiseur .trk capte la molette (scroll) et Esc/T (fermeture).
      g_trkScroll -= GetMouseWheelMove() * 40.0f;
      if (typed['t'] || IsKeyPressed(KEY_ESCAPE))
        showTrk = false;
    } else {
      if (overView) {
        // Zoom molette centre sur le curseur
        float wheel = GetMouseWheelMove();
        if (wheel != 0) {
          Vector2 worldBefore = GetScreenToWorld2D(mouse, cam);
          cam.zoom *= (wheel > 0) ? 1.1f : (1.0f / 1.1f);
          if (cam.zoom < 1e-6f)
            cam.zoom = 1e-6f;
          Vector2 worldAfter = GetScreenToWorld2D(mouse, cam);
          cam.target.x += worldBefore.x - worldAfter.x;
          cam.target.y += worldBefore.y - worldAfter.y;
        }
        // Pan : bouton droit ou milieu
        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT) || IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
          Vector2 d = GetMouseDelta();
          cam.target.x -= d.x / cam.zoom;
          cam.target.y -= d.y / cam.zoom;
        }
        // Selection : clic gauche
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
          pickGrain(GetScreenToWorld2D(mouse, cam));
        }
      }

      // Raccourcis clavier. Les lettres passent par typed[] (layout-aware) pour
      // coller au clavier physique (AZERTY) ; les fleches restent en scancode.
      if (IsKeyPressed(KEY_RIGHT))
        gotoDIC(g_dicNum + (iinc != 0 ? iinc : 1));
      if (IsKeyPressed(KEY_LEFT))
        gotoDIC(g_dicNum - (iinc != 0 ? iinc : 1));
      if (typed['f'])
        fitView(cam, sw, sh);
      if (typed['p'])
        showPattern = !showPattern;
      if (typed['z'])
        showSearchZones = !showSearchZones;
      if (typed['g'])
        showGrains = !showGrains;
      if (typed['r'])
        showTrajectory = !showTrajectory;
      if (typed['t'])
        showTrk = true;
      if (typed['b']) {
        showBackground = !showBackground;
        if (showBackground)
          loadBackground(g_dicNum);
        else
          unloadBackground();
      }
      if (IsKeyPressed(KEY_ESCAPE))
        break; // ESC quitte l'application
    }

    // ---- Dessin ----
    BeginDrawing();
    ClearBackground(COL_BG);

    BeginScissorMode(0, 0, viewW, sh);
    BeginMode2D(cam);

    if (showBackground && bgLoaded) {
      // La thumbnail couvre toute l'image : echelle = image_div pour revenir aux
      // coordonnees pixel pleine resolution dans lesquelles sont exprimes les grains.
      DrawTextureEx(bgTexture, {0, 0}, 0.0f, (float)bgLoadedDiv, WHITE);
    }
    drawGrains(cam.zoom);
    if (showTrajectory)
      drawTrajectories(cam.zoom);
    if (showPattern)
      drawPatterns();
    if (showSearchZones)
      drawSearchZones(cam.zoom);

    EndMode2D();
    EndScissorMode();

    // Colorbar du champ courant, en haut a droite de la vue (rien si "grey")
    drawColorBar(viewW);

    // ---- Bandeaux d'info flottants (cartes arrondies) ----
    {
      const char *info = g_fromDicFile ? TextFormat("DIC %d    |    %d grains", g_dicNum, (int)g_grains.size())
                                       : TextFormat("image %d    |    %d grains", g_dicNum, (int)g_grains.size());
      Vector2 ts = MeasureTextEx(uiFont, info, 20, 0);
      DrawRectangleRounded(Rectangle{10, 10, ts.x + 28, 38}, 0.35f, 8, Fade(COL_PANEL, 0.94f));
      uiText(info, 24, 19, 20, COL_TEXT);
    }

    // Panneau detaille du grain selectionne (toutes les infos disponibles)
    drawInspector();
    if (showSearchZones)
      drawSearchZonesLegend(sh);

    // ================= Panneau de controle (raygui) =================
    DrawRectangle(viewW, 0, PANEL_W, sh, COL_PANEL);
    DrawLine(viewW, 0, viewW, sh, COL_PANEL_EDGE); // separateur vertical

    const float px = (float)(viewW + 18);
    const float pw = (float)(PANEL_W - 36);
    const float BTN_H = (float)cfgBtnH; // hauteur des boutons / toggles / combos (configurable)
    const float SLD_H = (float)cfgSldH; // hauteur des sliders (configurable)
    const float LBL_H = (float)cfgLblH; // avance apres une ligne de texte (configurable)
    float y = 16;

    // En-tete
    uiText("D-TRACKER", px, y, 23, COL_TEXT);
    uiText("GUI", px + MeasureTextEx(uiFont, "D-TRACKER ", 23, 0).x, y + 5, 15, COL_ACCENT);
    y += 31;
    DrawLine((int)px, (int)y, (int)(px + pw), (int)y, COL_PANEL_EDGE);
    y += 3;

    // Helper de section : titre accentue + filet separateur
    auto section = [&](const char *title) {
      y += 8;
      uiText(title, px, y, 13, COL_ACCENT);
      y += 15;
      DrawLine((int)px, (int)y, (int)(px + pw), (int)y, COL_PANEL_EDGE);
      y += 6;
    };

    // ---- Navigation ----
    section("NAVIGATION");
    if (GuiButton(Rectangle{px, y, pw / 2 - 5, BTN_H}, "#118#  Prev"))
      gotoDIC(g_dicNum - (iinc != 0 ? iinc : 1));
    if (GuiButton(Rectangle{px + pw / 2 + 5, y, pw / 2 - 5, BTN_H}, "#119#  Next"))
      gotoDIC(g_dicNum + (iinc != 0 ? iinc : 1));
    y += BTN_H + 4;
    if (GuiButton(Rectangle{px, y, pw, BTN_H}, "#106#  Fit view  (F)"))
      fitView(cam, sw, sh);
    y += BTN_H + 4;

    // ---- Coloring ----
    section("COLORING");
    // Champ de coloration : grille de boutons (toutes les options visibles, 1 clic),
    // construite une fois a partir de COLOR_FIELDS (2 colonnes, ';' = colonne, '\n' = rangee).
    static const std::string fieldGrid = [] {
      std::string s;
      for (int i = 0; i < COLOR_FIELD_COUNT; ++i) {
        if (i > 0)
          s += (i % 2 == 0) ? '\n' : ';';
        s += COLOR_FIELDS[i].name;
      }
      return s;
    }();
    const int FGPAD = 2; // = GROUP_PADDING force dans applyCleanTheme
    const float FGH = 21;
    int fieldRows = (COLOR_FIELD_COUNT + 1) / 2;
    GuiToggleGroup(Rectangle{px, y, (pw - FGPAD) / 2, FGH}, fieldGrid.c_str(), &colorMode);
    y += fieldRows * FGH + (fieldRows - 1) * FGPAD + 6;

    // Borne(s) de la table de couleurs, selon le champ choisi.
    const ColorField &cf = COLOR_FIELDS[colorMode];
    if (cf.ncc) {
      uiText(TextFormat("NCC min  %.2f", colNCCMin), px, y, 14, COL_MUTED);
      y += LBL_H;
      GuiSliderBar(Rectangle{px, y, pw, SLD_H}, NULL, NULL, &colNCCMin, 0.0f, 0.975f);
      y += LBL_H;
    } else if (cf.member) {
      float lo, hi;
      colorRange(cf, lo, hi);
      uiText(TextFormat("range  [%.3g, %.3g]  auto", lo, hi), px, y, 14, COL_MUTED);
      y += LBL_H;
    }
    GuiToggle(Rectangle{px, y, pw, BTN_H}, showGrains ? "#44#  Hide grains  (G)" : "#44#  Show grains  (G)",
              &showGrains);
    y += BTN_H + 4;
    if (showGrains) {
      // L'opacite ne concerne que le remplissage des disques (pas la croix).
      uiText(TextFormat("Fill opacity  %.0f%%", grainAlpha * 100.0f), px, y, 14, COL_MUTED);
      y += LBL_H;
      GuiSliderBar(Rectangle{px, y, pw, SLD_H}, NULL, NULL, &grainAlpha, 0.0f, 1.0f);
      y += LBL_H;
    }
    GuiToggle(Rectangle{px, y, pw, BTN_H}, showPattern ? "#97#  Hide pattern  (P)" : "#97#  Show pattern  (P)",
              &showPattern);
    y += BTN_H + 4;
    GuiToggle(Rectangle{px, y, pw, BTN_H},
              showSearchZones ? "#98#  Hide search zones  (Z)" : "#98#  Show search zones  (Z)", &showSearchZones);
    y += BTN_H + 4;

    // ---- Positions (ref / precedente / actuelle) ----
    section("POSITIONS");
    GuiToggle(Rectangle{px, y, pw, BTN_H},
              showTrajectory ? "#90#  Hide ref/prev  (R)" : "#90#  Show ref/prev  (R)", &showTrajectory);
    y += BTN_H + 4;
    if (!showTrajectory)
      GuiSetState(STATE_DISABLED); // sous-option inactive tant que rien n'est affiche
    GuiToggle(Rectangle{px, y, pw, BTN_H}, "Selected grain only", &trajSelectedOnly);
    GuiSetState(STATE_NORMAL);
    y += BTN_H + 4;
    // Legende : pastille coloree + libelle pour chaque position
    float lx = px;
    auto legendDot = [&](Color col, const char *lab) {
      DrawCircleV(Vector2{lx + 5, y + 8}, 5, col);
      uiText(lab, lx + 14, y, 14, COL_MUTED);
      lx += 14 + MeasureTextEx(uiFont, lab, 14, 0).x + 14;
    };
    legendDot(COL_POS_REF, "ref");
    legendDot(COL_POS_PREV, "prev");
    legendDot(COL_POS_CUR, "cur");
    y += LBL_H;

    // ---- Background image ----
    section("BACKGROUND IMAGE");
    bool prevBg = showBackground;
    GuiToggle(Rectangle{px, y, pw, BTN_H}, showBackground ? "#12#  Hide image  (B)" : "#12#  Show image  (B)",
              &showBackground);
    if (showBackground != prevBg) {
      if (showBackground)
        loadBackground(g_dicNum);
      else
        unloadBackground();
    }
    y += BTN_H + 4;
    uiText(TextFormat("Resolution divider  %d", gui_image_div), px, y, 14, COL_MUTED);
    y += LBL_H;
    float fdiv = (float)gui_image_div;
    GuiSliderBar(Rectangle{px, y, pw, SLD_H}, NULL, NULL, &fdiv, 1.0f, 32.0f);
    int newDiv = (int)(fdiv + 0.5f);
    if (newDiv != gui_image_div) {
      gui_image_div = newDiv;
      if (showBackground)
        loadBackground(g_dicNum); // recharge a la nouvelle resolution
    }
    y += LBL_H;

    // ---- Computation ----
    section("COMPUTATION");
    uiText(TextFormat("procedure: %s", procedure.c_str()), px, y, 14, COL_MUTED);
    y += 20;
    if (GuiButton(Rectangle{px, y, pw, 26}, TextFormat("#131#  %s", procedureLabel())))
      runCorrelation();
    y += 30;
    if (GuiButton(Rectangle{px, y, pw, BTN_H}, showTrk ? "#10#  Hide .trk  (T)" : "#10#  View .trk  (T)"))
      showTrk = !showTrk;
    y += BTN_H + 4;

    // ---- Settings ----
    section("SETTINGS");
    if (GuiButton(Rectangle{px, y, pw, BTN_H}, "#2#  Save settings  (.ini)")) {
      cfgWinW = sw; // capture la taille courante de la fenetre
      cfgWinH = sh;
      saveConfig();
      g_resultLines = {TextFormat("Settings saved -> %s", CONFIG_FILE)};
      g_resultUntil = GetTime() + 2.5;
    }
    y += BTN_H + 4;

    // Bandeau recapitulatif du dernier Run (temporaire)
    drawResultBanner(viewW, sh);

    // Visualiseur .trk : overlay au premier plan (par-dessus tout le reste)
    drawTrkViewer(viewW, sh);

    EndDrawing();
  }

  UnloadFont(uiFont);

  unloadBackground();
  CloseWindow();
  return 0;
}
