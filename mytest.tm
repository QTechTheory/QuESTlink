
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
:Evaluate: AnotherTest::usage = "AnotherTest[] gives summin."


:Begin:
:Function:       returnQureg
:Pattern:        ReturnQureg[]
:Arguments:      {}
:ArgumentTypes:  {}
:ReturnType:     Manual
:End:
:Evaluate: ReturnQureg::usage = "ReturnQureg[] gives summin."


:Begin:
:Function:       giveQureg
:Pattern:        GiveQureg[q_Association]
:Arguments:      { q }
:ArgumentTypes:  { Manual }
:ReturnType:     Manual
:End:
:Evaluate: GiveQureg::usage = "GiveQureg[list] gives summin."