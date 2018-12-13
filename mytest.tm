
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