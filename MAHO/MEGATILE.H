// Andriy_Zatserklyaniy@fnal.gov

struct Point;
struct Side;
struct Groove;

template <class T> class Queue;
class TileData;

// for graphics
class TileSide;

enum EndStatus {tofly, tocut, tosemihole, tohole, torotate};
enum EndCutter {noncut, cutcut, cylcut, balcut};

// general-purpose functions

struct End;
struct EndGroove;

int ToCut(const Queue<End> *list);
int ToCut(const Queue<EndGroove> *list);

// struct End

struct End {
private:
  static short EndAnchor;
public:
  static void* Key() {return (void*) &EndAnchor;}
public:
  double X,Y;
  double Dcutter;
  EndCutter cutype;
  EndStatus status;
  TileSide *tile;

  End() {X=Y=Dcutter=0.; cutype=cutcut; status=tofly; tile=NULL;}
  End(const End& end) {
    X=end.X; Y=end.Y;
    Dcutter=end.Dcutter;
    cutype=end.cutype;
    status=end.status;
    tile = end.tile;
  }
  End(const Point& p);
  End(const Point& p, double dcutter);
  End(const Point& p, EndStatus astatus);
  End(const Point& p, double dcutter, EndStatus stat);
  End(const Point& p, double dcutter, EndStatus stat, TileSide *master);
  End(const Point& p, double dcutter, EndCutter ecut, EndStatus stat);
  End(const Point& p, double dcutter, EndCutter ecut, EndStatus stat,
      TileSide *master);
  virtual void* GetKey() {return (void*) &EndAnchor;}
  virtual void Show() const;
};

struct EndGroove: public End {
private:
  static short EndGrooveAnchor;
public:
  static void* Key() {return (void*) &EndGrooveAnchor;}
public:
  Queue<End> *P;

  virtual void* GetKey() {return (void*) &EndGrooveAnchor;}
  virtual void Show() const;
  EndGroove() : End() {P = new Queue<End>;}
 ~EndGroove() {delete P; P=NULL;}
};

// struct is used to select tile to cut
struct R2TileSidePoint;

class Megatile {
  // static variables
public:
  static double Dcut;
  static double Dhole;
  static char   task_file[33];
  static char   path_file[33];
  static char   groove_file[33];
  static double SheetX1;
  static double SheetY1;
  static double SheetX2;
  static double SheetY2;
  // shift with respect to (SheetX1,SheetY1)
  static double xShift;
  static double yShift;

  static int TrueAngle;

  static int SqueezeParts;
  static int BackSide;
  static int down_cut;
  static int down_todo;
  // extra space to make down cut
  static double ddown;

  // safety gates
  static Queue<Side> *gates;

  // tiles
  static Queue<char> *NameList;

  // start after break
  static char continue_tile[33];
  static int  continue_level;

  // optional show of turned contour
  static int ShowTurnedContour;
public:
  static void WhoAlive();
  static void WhoAlive(char *mess);

protected:
  // center of rotation
  const Point O;

private:
  // Lists of pointers to tiles in the right and left parts of sheet
  // which are use by Squeeze() and old MakeDownCut()
  Queue<TileData> *RightPart;
  Queue<TileData> *LeftPart;

  Side *continue_cut;
private:
  void MakeName(char* name, int megatile, int megaside, int number);
  int ExistName(char *name, Queue<TileData> *tds);
  int IndexName(char *name, Queue<TileData> *tds);
  TileData* tdName(char *name, Queue<TileData> *tds);
  int SelectTileToCut(R2TileSidePoint& tsp);
  int CheckStatus(TileData *tile);
  Side* DownSide(Queue<TileData> *td);
  void MakeDownCut();

  Point DownRightPoint;
  Point    UpLeftPoint;

  Queue<Side> *cuts;
  Queue<Side> *down_cuts;
protected:
  Queue<TileData> *tds;

  Queue<End> *ends;
  Queue<End> *down_ends;
  Queue<EndGroove> *groove_ends;
protected:
  int nopath;
public:
  Megatile(double theOx, double theOy);
  ~Megatile();
  void Insert(char *tname, double x, double y, double thdeg, int inv);
  void Insert(char *tname, double thdeg, int inv);
  void Insert(char *tname, int inv);
  void Insert(char *tname);

public:
  //void MakePath();
protected:
  void MakePath();
  void Sort();
  double dxRightPart(const Point& point);
  double dxLeftPart(const Point& point);
  double dxPart(const Queue<TileData> *list, const Point& point);
  double dxPart(const Queue<TileData> *list, const Point& point, double x0);
  TileData* ClosestTileX (const Queue<TileData> *list, const Point& corner);
  Point ClosestX (const Queue<TileData> *list, const Point& corner);
  Point ClosestY (const Queue<TileData> *list, const Point& corner);
  Point ClosestXY(const Queue<TileData> *list, const Point& corner);
  Point ClosestYX(const Queue<TileData> *list, const Point& corner);
  Point ContourClosestX (const Queue<TileData> *list, const Point& corner);
  void Squeeze();
  void DefineLimitPoints();
  void Arrange();
  void ArrangeMegaSide1();
  void ArrangeMegaSide1Right();
  void ArrangeMegaSide1Left();
  void ArrangeMegaSide2();
  void ArrangeMegaSide2Right();
  void ArrangeMegaSide2Left();
  int  SheetIncludes(TileSide *tile);
  int  InSheetWorkArea();
  int  TilesOverlap();
  int  Fit();


  void CutPath();
  void DownEnds();
  TileSide* Master(Side *side);
  void AddCutToEnds(Side *cut, Queue<End> *endlist, TileSide *master);
  void AddCutToEnds(Side *cut, Queue<End> *endlist) {
    TileSide *master = Master(cut);
    AddCutToEnds(cut,endlist,master);
  }
  void CuttingEnds();
  void GroovingEnds();
  void AddGrooves(Queue<TileData> *list, EndCutter ecut);

  void PrintCuts();
  void PrintEnds();
  void WriteEnds();
  void WriteGrooves();
  void PrintCutsEnds();

  void Rotate(const Point& O, double theta);
  void TurnBack();

  // for graphics package
  char* GetName(TileSide *tile);
  Point TileUpLeftPoint(TileSide *tile);
  /*
  int TilesGetCount();
  int SidesGetCount(int ntile);
  Queue<Side>*  GetContour(int ntile);
  Queue<Side>*  GetTurnedContour(int ntile);
  Queue<Side>*  GetSides(int ntile);
  Groove* GetGroove(int ntile);
  char* GetName(int ntile);
  Point UpLeftPoint(int ntile);
  */
public:
  void Area(double& xmin, double& xmax, double& ymin, double& ymax);
  void Show() const;
};
