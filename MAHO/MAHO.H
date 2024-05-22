// Andriy_Zatserklyaniy@fnal.gov

template <class T> class Queue;
class Megatile;

// cutter current status
enum CutterStatusType {cut,fly};
enum G0G1statusType   {G00,G0,G1};

class MAHO: public Megatile {
public:
  // static variables
  static int G28on;

  //Machine work area
  static double MAHO_XMIN;
  static double MAHO_XMAX;
  static double MAHO_ZMIN;
  static double MAHO_ZMAX;

  // WorkPiece zero
  static double MAHO_X0;
  static double MAHO_Y0;
  static double MAHO_Z0;
  // The center of rotation
  static double MAHO_OX;
  static double MAHO_OZ;

  static double ChangeCutter_Y;
  static double ChangeCutter_Z;

  static int NumberDcut;
  static int NumberDsemihole;
  static int NumberDhole;
  static int NumberDcyl;
  static int NumberDbal;
  static char code_file[33];
  static double hfly;
  // the height of cutter top position
  static double htop;
  // rotation
  static long rotation;
  // direction of rotation
  static int dir_rotation;

  static Queue<char>   *comment;

  // sheet thickness
  static double         thickness;
  // the min cutter height over the table
  static double         minheight;
  // correct the last cut depth for minheight-thickness
  static int            correct_depth;

  static Queue<double> *depth;
  static Queue<long>   *feed;
  static Queue<double> *cdepth;
  static Queue<long>   *cfeed;
  static Queue<double> *bdepth;
  static Queue<long>   *bfeed;
  // No. of prigram ( >= 9000 )
  static long ProgramNumber;
  // No. of prigram
  static long SubNumber;
  static long LineInc;

private:
  // cutter_status == fly corresponds to height hfly,
  //                  cut corresponds to current depth
  CutterStatusType CutterStatus;

  G0G1statusType G0G1status;
  EndCutter CurrentCutter;
  double    CurrentDepth;
  long      CurrentFeed;
  double    CurrentAngle;
  TileSide *CurrentTile;

  int NoRotation;
  int NoFeed;

  // string number
  int nmain;
  int nsub;

  // buffer
  char BUF[255];

  Queue<char> **part, *code, *main, *sub;

  int CurrentLine() {return (part==&main)? nmain-LineInc: nsub-LineInc;}
  int NextLine() {return (part==&main)? nmain: nsub;}

  // put string to current part
  void PUT();

  void MainHeader();
  void MakeMain();

  void MakeSubMain();
  void SubHeader();
  void MakeSub();
  // MAHO commands
  void Stop();
  void Pause(char *mess);
  void Rotate(const double angle);
  void LiftCutterFast();
  // commands which use G1
  void LiftCutter();
  void DownCutter();
  void MoveHeight(double h);
  void Cut(const End *end);

  void SetCutterAtEnd(const End *end);
  void MakeHole(const End *end);
  void Fly(const End *end);
  void CutGroove(const EndGroove *groove);

  void SetDepth(double thedepth)
  {
    if (CurrentDepth != thedepth) {
      CurrentDepth = thedepth;
    }
  }
  void SetFeed(long thefeed)
  {
    if (CurrentFeed != thefeed) {
      CurrentFeed=thefeed;
      NoFeed=1;
    }
  }

  void CutCache(const Queue<End> *path);
  void CutPath(const Queue<End> *path);
  void CutGrooves();

public:
  MAHO();
  ~MAHO();
  void Insert(char *name, double x, double y, double thdeg, int inv) {
    Megatile::Insert(name,x,y,thdeg,inv);
  }
  void Insert(char *tname, int inv) {
    Megatile::Insert(tname,inv);
  }
  void Insert(char *tname) {
    Megatile::Insert(tname);
  }
  void MakePath() {
    Megatile::MakePath();
    Megatile::WriteEnds();
    //Megatile::PrintCutsEnds();
  }
  void MakeCode();
  void PrintCode();
  void WriteCode();
  //void Graph() {Megatile::Graph();}
};
