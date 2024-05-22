// Andriy_Zatserklyaniy@fnal.gov

struct viewporttype;
struct textsettingstype;

// general-purpose function
void DashedArc(int x, int y, int theta1, int theta2, int R, int dtheta);

template <class T> class Queue;

//class Megatile;
class TileData;
class MAHO;

struct Window {
private:
public:
  /*   Global variables   */
  static struct viewporttype VPORTCURR;
  static struct textsettingstype TEXTSETCURR;
  static int COLORCURR;
public:
  static int hText();
  static int wText();

protected:
  int X1,Y1,X2,Y2,Clip;
  int dX,dY;
public:
  Window(int x1, int y1, int x2, int y2, int clip=1);
  Window();
  void DefineWindow(int x1, int y1, int x2, int y2, int clip=1);
  void SaveCurrentWindow();
  void RestoreSavedWindow();
  int GetX1() {return X1;}
  int GetX2() {return X2;}
  int GetY1() {return Y1;}
  int GetY2() {return Y2;}
  virtual void Set();
  virtual void Clear();
  virtual void Show();
};

struct MenuWindow: public Window {
  int xcurr, ycurr;
  int TextColor;
  Queue<char> *items;

   MenuWindow(int x1, int y1, int x2, int y2, int clip, int color);
   MenuWindow();
  ~MenuWindow();
  void DefineMenuWindow(int x1, int y1, int x2, int y2, int clip, int color);
  void Show();
  void Set();
  void writeln(char *item);
  void Add(char *Item);
};

//class GManager: public Megatile {
class GManager: public MAHO {
public:
  static int graphics;

  static long CuttingDelay;  // in milliseconds
  static long FlyingDelay;   // in milliseconds
  static long GroovingDelay; // in milliseconds

  static int Color;
  static int BkColor;
  static int TextColor;
  static int MessageColor;
  static int CutColor;
  static int FlyColor;
  static int ContourColor;
  static int SideColor;
  static int CylGrooveColor;
  static int BallGrooveColor;
  static int GatesColor;
  static int SheetColor;
  static int WorkAreaColor;
  static int FinishColor;
  static int Rfinish;
  static int FinishSnow;

  int MAXX, MAXY, MAXCOLOR;

  // global used by cut functions
  End *FROMend;

  Window     *work;
  MenuWindow *main;
  Window     *mess;
  Window     *entire;

  int code_written;
  int Ruler;
  int MachineArea;
  int SheetArea;
  int Names;
  int hitCR;

   GManager();
  ~GManager();
  void InitGraph();
  void Area();
  void GraphName(TileSide *tile);
  void GraphNames();
  void GraphContour(TileSide *tile);
  void GraphTurnedContour(TileSide *tile);
  void GraphSides(TileSide *tile);
  int  GraphGroove(Groove* groove, int nparts);
  int  GraphCylGroove(TileSide *tile);
  int  GraphBalGroove(TileSide *tile);
  void GraphGrooveContour(TileSide *tile, int color);
  void GraphTile(TileSide *tile);
  void GraphTiles();
  void GraphGates();
  void GraphRuler();
  void GraphSheet();
  int  GraphCutting();
  void GraphCut(End *end);
  void GraphFly(End *end);
  void GraphMakeHole(End *end);
  void GraphLiftCutter();
  void Heap();
  void Graph();
  void WriteCode();
  int Xscr(double x);
  int Yscr(double y);
  void Xtick(double x, double y, int h);
  void Xlabl(double x, double y);
  void Ytick(double x, double y, int h);
  void Ylabl(double x, double y);
  void Slower();
  void Faster();
  int  Pause();
  void Quit();
  void Snow();
  static int Graphics() {return graphics;}
  static void CloseGraph();

  void DefineWindows();
  void Message(char *str);
  void Message(char *str1, char *str2);
  void Menu();

  void Insert(char *name, double x, double y, double thdeg, int inv) {
    Megatile::Insert(name,x,y,thdeg,inv);
  }
  void Insert(char *tname, int inv) {
    Megatile::Insert(tname,inv);
  }

  void Insert(char *tname, double thdeg, int inv) {
    Megatile::Insert(tname,thdeg,inv);
  }

  void Insert(char *tname) {
    Megatile::Insert(tname);
  }

  void MakePath() {
    Megatile::MakePath();
    Megatile::WriteEnds();
    //Megatile::PrintCutsEnds();
  }
};

