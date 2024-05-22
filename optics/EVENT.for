*                                         @METAGS event
*                                         09-23-98 07:58pm
*--------------- event ---------------
*
      function event(dummy)
      include 'include.inc'
      print*, 'CFILE =', CFILE
      print*, 'IDNEVT =', IDNEVT
      ncur = IDNEVT
      IDNEVT = 1
      print*, 'IDNEVT =', IDNEVT
      print*, 'Rbox =', Rbox
      IDNEVT = ncur
      print*, 'IDNEVT =', IDNEVT
      event = 1.
      END
