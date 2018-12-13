
:Begin:
:Function:       mytest
:Pattern:        MyTest[i_Integer, j_Integer]
:Arguments:      { i, j }
:ArgumentTypes:  { Integer, Integer }
:ReturnType:     Integer
:End:
:Evaluate: MyTest::usage = "MyTest[x, y] gives the sum of two machine integers x and y."


:Begin:
:Function:       anothertest
:Pattern:        AnotherTest[]
:Arguments:      {}
:ArgumentTypes:  {}
:ReturnType:     Real
:End:
:Evaluate: MyTest::usage = "AnotherTest[] gives summin."


:Begin:
:Function:       getQureg
:Pattern:        GetQureg[]
:Arguments:      {}
:ArgumentTypes:  {}
:ReturnType:     Manual
:End:
:Evaluate: MyTest::usage = "GetQureg[] gives summin."


:Begin:
:Function:       giveQureg
:Pattern:        GiveQureg[i_Integer, j_Integer]
:Arguments:      { i, j }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: MyTest::usage = "GiveQureg[] gives summin."