/*
 * This file automatically produced by ./WSTPlibs/wsprep from:
 *	quest_templates.tm
 * mprep Revision 18 Copyright (c) Wolfram Research, Inc. 1990-2013
 */

#define MPREP_REVISION 18
#include "wstp.h"


int PREPAbort = 0;
int PREPDone  = 0;
long PREPSpecialCharacter = '\0';

WSLINK stdlink = 0;
WSEnvironment stdenv = 0;
WSYieldFunctionObject stdyielder = (WSYieldFunctionObject)0;
WSMessageHandlerObject stdhandler = (WSMessageHandlerObject)0;


int wrapper_createQureg P(( int _tp1));

static int _tr0( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSNewPacket(wslp) ) goto L1;

	_rp0 = wrapper_createQureg(_tp1);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L1: 
L0:	return res;
} /* _tr0 */


int wrapper_createDensityQureg P(( int _tp1));

static int _tr1( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSNewPacket(wslp) ) goto L1;

	_rp0 = wrapper_createDensityQureg(_tp1);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L1: 
L0:	return res;
} /* _tr1 */


int wrapper_destroyQureg P(( int _tp1));

static int _tr2( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSNewPacket(wslp) ) goto L1;

	_rp0 = wrapper_destroyQureg(_tp1);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L1: 
L0:	return res;
} /* _tr2 */


int wrapper_initZeroState P(( int _tp1));

static int _tr3( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSNewPacket(wslp) ) goto L1;

	_rp0 = wrapper_initZeroState(_tp1);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L1: 
L0:	return res;
} /* _tr3 */


int wrapper_initPlusState P(( int _tp1));

static int _tr4( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSNewPacket(wslp) ) goto L1;

	_rp0 = wrapper_initPlusState(_tp1);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L1: 
L0:	return res;
} /* _tr4 */


int wrapper_initClassicalState P(( int _tp1, int _tp2));

static int _tr5( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSGetInteger( wslp, &_tp2) ) goto L1;
	if ( ! WSNewPacket(wslp) ) goto L2;

	_rp0 = wrapper_initClassicalState(_tp1, _tp2);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L2: L1: 
L0:	return res;
} /* _tr5 */


int wrapper_initPureState P(( int _tp1, int _tp2));

static int _tr6( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSGetInteger( wslp, &_tp2) ) goto L1;
	if ( ! WSNewPacket(wslp) ) goto L2;

	_rp0 = wrapper_initPureState(_tp1, _tp2);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L2: L1: 
L0:	return res;
} /* _tr6 */


int wrapper_initStateFromAmps P(( int _tp1, double * _tp2, long _tpl2, double * _tp3, long _tpl3));

static int _tr7( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	double * _tp2;
	long _tpl2;
	double * _tp3;
	long _tpl3;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSGetRealList( wslp, &_tp2, &_tpl2) ) goto L1;
	if ( ! WSGetRealList( wslp, &_tp3, &_tpl3) ) goto L2;
	if ( ! WSNewPacket(wslp) ) goto L3;

	_rp0 = wrapper_initStateFromAmps(_tp1, _tp2, _tpl2, _tp3, _tpl3);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L3:	WSReleaseReal64List(wslp, _tp3, _tpl3);
L2:	WSReleaseReal64List(wslp, _tp2, _tpl2);
L1: 
L0:	return res;
} /* _tr7 */


int wrapper_cloneQureg P(( int _tp1, int _tp2));

static int _tr8( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSGetInteger( wslp, &_tp2) ) goto L1;
	if ( ! WSNewPacket(wslp) ) goto L2;

	_rp0 = wrapper_cloneQureg(_tp1, _tp2);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L2: L1: 
L0:	return res;
} /* _tr8 */


int wrapper_applyOneQubitDepolariseError P(( int _tp1, int _tp2, double _tp3));

static int _tr9( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _tp2;
	double _tp3;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSGetInteger( wslp, &_tp2) ) goto L1;
	if ( ! WSGetReal( wslp, &_tp3) ) goto L2;
	if ( ! WSNewPacket(wslp) ) goto L3;

	_rp0 = wrapper_applyOneQubitDepolariseError(_tp1, _tp2, _tp3);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L3: L2: L1: 
L0:	return res;
} /* _tr9 */


int wrapper_applyTwoQubitDepolariseError P(( int _tp1, int _tp2, int _tp3, double _tp4));

static int _tr10( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	double _tp4;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSGetInteger( wslp, &_tp2) ) goto L1;
	if ( ! WSGetInteger( wslp, &_tp3) ) goto L2;
	if ( ! WSGetReal( wslp, &_tp4) ) goto L3;
	if ( ! WSNewPacket(wslp) ) goto L4;

	_rp0 = wrapper_applyTwoQubitDepolariseError(_tp1, _tp2, _tp3, _tp4);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L4: L3: L2: L1: 
L0:	return res;
} /* _tr10 */


int wrapper_applyOneQubitDephaseError P(( int _tp1, int _tp2, double _tp3));

static int _tr11( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _tp2;
	double _tp3;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSGetInteger( wslp, &_tp2) ) goto L1;
	if ( ! WSGetReal( wslp, &_tp3) ) goto L2;
	if ( ! WSNewPacket(wslp) ) goto L3;

	_rp0 = wrapper_applyOneQubitDephaseError(_tp1, _tp2, _tp3);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L3: L2: L1: 
L0:	return res;
} /* _tr11 */


int wrapper_applyTwoQubitDephaseError P(( int _tp1, int _tp2, int _tp3, double _tp4));

static int _tr12( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	double _tp4;
	int _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSGetInteger( wslp, &_tp2) ) goto L1;
	if ( ! WSGetInteger( wslp, &_tp3) ) goto L2;
	if ( ! WSGetReal( wslp, &_tp4) ) goto L3;
	if ( ! WSNewPacket(wslp) ) goto L4;

	_rp0 = wrapper_applyTwoQubitDephaseError(_tp1, _tp2, _tp3, _tp4);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);
L4: L3: L2: L1: 
L0:	return res;
} /* _tr12 */


double wrapper_calcProbOfOutcome P(( int _tp1, int _tp2, int _tp3));

static int _tr13( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	double _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSGetInteger( wslp, &_tp2) ) goto L1;
	if ( ! WSGetInteger( wslp, &_tp3) ) goto L2;
	if ( ! WSNewPacket(wslp) ) goto L3;

	_rp0 = wrapper_calcProbOfOutcome(_tp1, _tp2, _tp3);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutReal( wslp, _rp0);
L3: L2: L1: 
L0:	return res;
} /* _tr13 */


double wrapper_calcFidelity P(( int _tp1, int _tp2));

static int _tr14( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	int _tp2;
	double _rp0;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSGetInteger( wslp, &_tp2) ) goto L1;
	if ( ! WSNewPacket(wslp) ) goto L2;

	_rp0 = wrapper_calcFidelity(_tp1, _tp2);

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutReal( wslp, _rp0);
L2: L1: 
L0:	return res;
} /* _tr14 */


void internal_applyCircuit P(( int _tp1));

static int _tr15( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;

	internal_applyCircuit(_tp1);

	res = 1;
 
L0:	return res;
} /* _tr15 */


void internal_getStateVec P(( int _tp1));

static int _tr16( WSLINK wslp)
{
	int	res = 0;
	int _tp1;
	if ( ! WSGetInteger( wslp, &_tp1) ) goto L0;
	if ( ! WSNewPacket(wslp) ) goto L1;

	internal_getStateVec(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr16 */


int callable_destroyAllQuregs P(( void));

static int _tr17( WSLINK wslp)
{
	int	res = 0;
	int _rp0;
	if ( ! WSNewPacket(wslp) ) goto L0;

	_rp0 = callable_destroyAllQuregs();

	res = PREPAbort ?
		WSPutFunction( wslp, "Abort", 0) : WSPutInteger( wslp, _rp0);

L0:	return res;
} /* _tr17 */


void callable_getAllQuregs P(( void));

static int _tr18( WSLINK wslp)
{
	int	res = 0;
	if ( ! WSNewPacket(wslp) ) goto L0;
	if( !wslp) return res; /* avoid unused parameter warning */

	callable_getAllQuregs();

	res = 1;

L0:	return res;
} /* _tr18 */


static struct func {
	int   f_nargs;
	int   manual;
	int   (*f_func)P((WSLINK));
	const char  *f_name;
	} _tramps[19] = {
		{ 1, 0, _tr0, "wrapper_createQureg" },
		{ 1, 0, _tr1, "wrapper_createDensityQureg" },
		{ 1, 0, _tr2, "wrapper_destroyQureg" },
		{ 1, 0, _tr3, "wrapper_initZeroState" },
		{ 1, 0, _tr4, "wrapper_initPlusState" },
		{ 2, 0, _tr5, "wrapper_initClassicalState" },
		{ 2, 0, _tr6, "wrapper_initPureState" },
		{ 3, 0, _tr7, "wrapper_initStateFromAmps" },
		{ 2, 0, _tr8, "wrapper_cloneQureg" },
		{ 3, 0, _tr9, "wrapper_applyOneQubitDepolariseError" },
		{ 4, 0, _tr10, "wrapper_applyTwoQubitDepolariseError" },
		{ 3, 0, _tr11, "wrapper_applyOneQubitDephaseError" },
		{ 4, 0, _tr12, "wrapper_applyTwoQubitDephaseError" },
		{ 3, 0, _tr13, "wrapper_calcProbOfOutcome" },
		{ 2, 0, _tr14, "wrapper_calcFidelity" },
		{ 1, 2, _tr15, "internal_applyCircuit" },
		{ 1, 0, _tr16, "internal_getStateVec" },
		{ 0, 0, _tr17, "callable_destroyAllQuregs" },
		{ 0, 0, _tr18, "callable_getAllQuregs" }
		};

static const char* evalstrs[] = {
	"QuEST`CreateQureg::usage = \"CreateQureg[numQubits] returns the i",
	"d of a newly created remote statevector.\"",
	(const char*)0,
	"QuEST`CreateDensityQureg::usage = \"CreateDensityQureg[numQubits]",
	" returns the id of a newly created remote density matrix.\"",
	(const char*)0,
	"QuEST`Private`DestroyQuregInternal::usage = \"DestroyQuregInterna",
	"l[numQubits] frees the memory of the remote qureg associated wit",
	"h the given id.\"",
	(const char*)0,
	"QuEST`InitZeroState::usage = \"InitZeroState[qureg] returns a sta",
	"te in |0>.\"",
	(const char*)0,
	"QuEST`InitPlusState::usage = \"InitPlusState[qureg] returns a sta",
	"te in |+>.\"",
	(const char*)0,
	"QuEST`InitClassicalState::usage = \"InitClassicalState[qureg, ind",
	"] returns a state in basis state |ind>.\"",
	(const char*)0,
	"QuEST`InitPureState::usage = \"InitPureState[targetQureg, pureQur",
	"eg] puts targetQureg (statevec or density matrix) into the pureQ",
	"ureg (statevec) state.\"",
	(const char*)0,
	"QuEST`InitStateFromAmps::usage = \"InitStateFromAmps[qureg, reals",
	", imags] initialises the given qureg to have the supplied amplit",
	"udes.\"",
	(const char*)0,
	"QuEST`CloneQureg::usage = \"CloneQureg[dest, source] sets dest to",
	" be a copy of source.\"",
	(const char*)0,
	"QuEST`ApplyOneQubitDepolariseError::usage = \"ApplyOneQubitDepola",
	"riseError[qureg, qubit, prob] adds depolarising noise to density",
	" matrix qureg.\"",
	(const char*)0,
	"QuEST`ApplyTwoQubitDepolariseError::usage = \"ApplyTwoQubitDepola",
	"riseError[qureg, qb1, qb2 prob] adds depolarising noise to densi",
	"ty matrix qureg.\"",
	(const char*)0,
	"QuEST`ApplyOneQubitDephaseError::usage = \"ApplyOneQubitDephaseEr",
	"ror[qureg, qubit, prob] adds dephasing noise to density matrix q",
	"ureg.\"",
	(const char*)0,
	"QuEST`ApplyTwoQubitDephaseError::usage = \"ApplyTwoQubitDephaseEr",
	"ror[qureg, qb1, qb2 prob] adds dephasing noise to density matrix",
	" qureg.\"",
	(const char*)0,
	"QuEST`CalcProbOfOutcome::usage = \"CalcProbOfOutcome[qureg, qubit",
	", outcome] returns the probability of measuring qubit in the giv",
	"en outcome.\"",
	(const char*)0,
	"QuEST`CalcFidelity::usage = \"CalcFidelity[qureg1, qureg2] return",
	"s the fidelity between the given states.\"",
	(const char*)0,
	"QuEST`Private`ApplyCircuitInner::usage = \"ApplyCircuitInner[qure",
	"g, opcodes, ctrls, numCtrlsPerOptrls, targs, params] applies a c",
	"ircuit (decomposed into codes) to the given qureg.\"",
	(const char*)0,
	"QuEST`Private`GetStateVecInternal::usage = \"GetStateVecInternal[",
	"qureg] returns the underlying statevector associated with the gi",
	"ven qureg (flat, even for density matrices).\"",
	(const char*)0,
	"QuEST`DestroyAllQuregs::usage = \"DestroyAllQuregs[] destroys all",
	" remote quregs.\"",
	(const char*)0,
	"QuEST`GetAllQuregs::usage = \"GetAllQuregs[] returns all active q",
	"uregs.\"",
	(const char*)0,
	(const char*)0
};
#define CARDOF_EVALSTRS 19

static int _definepattern P(( WSLINK, char*, char*, int));

static int _doevalstr P(( WSLINK, int));

int  _PREPDoCallPacket P(( WSLINK, struct func[], int));


int WSInstall( WSLINK wslp)
{
	int _res;
	_res = WSConnect(wslp);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`CreateQureg[numQubits_Integer]", (char *)"{ numQubits }", 0);
	if (_res) _res = _doevalstr( wslp, 0);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`CreateDensityQureg[numQubits_Integer]", (char *)"{ numQubits }", 1);
	if (_res) _res = _doevalstr( wslp, 1);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`Private`DestroyQuregInternal[id_Integer]", (char *)"{ id }", 2);
	if (_res) _res = _doevalstr( wslp, 2);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`InitZeroState[qureg_Integer]", (char *)"{ qureg }", 3);
	if (_res) _res = _doevalstr( wslp, 3);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`InitPlusState[qureg_Integer]", (char *)"{ qureg }", 4);
	if (_res) _res = _doevalstr( wslp, 4);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`InitClassicalState[qureg_Integer, state_Integer]", (char *)"{ qureg, state }", 5);
	if (_res) _res = _doevalstr( wslp, 5);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`InitPureState[targetQureg_Integer, pureQureg_Integer]", (char *)"{ targetQureg, pureQureg }", 6);
	if (_res) _res = _doevalstr( wslp, 6);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`InitStateFromAmps[qureg_Integer, reals_List, imags_List]", (char *)"{ qureg, reals, imags }", 7);
	if (_res) _res = _doevalstr( wslp, 7);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`CloneQureg[target_Integer, source_Integer]", (char *)"{ target, source }", 8);
	if (_res) _res = _doevalstr( wslp, 8);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`ApplyOneQubitDepolariseError[qureg_Integer, qb_Integer, prob_Real]", (char *)"{ qureg, qb, prob }", 9);
	if (_res) _res = _doevalstr( wslp, 9);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`ApplyTwoQubitDepolariseError[qureg_Integer, qb1_Integer, qb2_Integer, prob_Real]", (char *)"{ qureg, qb1, qb2, prob }", 10);
	if (_res) _res = _doevalstr( wslp, 10);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`ApplyOneQubitDephaseError[qureg_Integer, qb_Integer, prob_Real]", (char *)"{ qureg, qb, prob }", 11);
	if (_res) _res = _doevalstr( wslp, 11);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`ApplyTwoQubitDephaseError[qureg_Integer, qb1_Integer, qb2_Integer, prob_Real]", (char *)"{ qureg, qb1, qb2, prob }", 12);
	if (_res) _res = _doevalstr( wslp, 12);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`CalcProbOfOutcome[qureg_Integer, qb_Integer, outcome_Integer]", (char *)"{ qureg, qb, outcome }", 13);
	if (_res) _res = _doevalstr( wslp, 13);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`CalcFidelity[qureg1_Integer, qureg2_Integer]", (char *)"{ qureg1, qureg2 }", 14);
	if (_res) _res = _doevalstr( wslp, 14);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`Private`ApplyCircuitInternal[qureg_Integer, opcodes_List, ctrls_List, numCtrlsPerOp_List, targs_List, params_List]", (char *)"{ qureg, opcodes, ctrls, numCtrlsPerOp, targs, params }", 15);
	if (_res) _res = _doevalstr( wslp, 15);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`Private`GetStateVecInternal[qureg_Integer]", (char *)"{ qureg }", 16);
	if (_res) _res = _doevalstr( wslp, 16);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`DestroyAllQuregs[]", (char *)"{ }", 17);
	if (_res) _res = _doevalstr( wslp, 17);
	if (_res) _res = _definepattern(wslp, (char *)"QuEST`GetAllQuregs[]", (char *)"{ }", 18);
	if (_res) _res = _doevalstr( wslp, 18);
	if (_res) _res = WSPutSymbol( wslp, "End");
	if (_res) _res = WSFlush( wslp);
	return _res;
} /* WSInstall */


int PREPDoCallPacket( WSLINK wslp)
{
	return _PREPDoCallPacket( wslp, _tramps, 19);
} /* PREPDoCallPacket */

/******************************* begin trailer ********************************/

#ifndef EVALSTRS_AS_BYTESTRINGS
#	define EVALSTRS_AS_BYTESTRINGS 1
#endif

#if CARDOF_EVALSTRS
static int  _doevalstr( WSLINK wslp, int n)
{
	long bytesleft, charsleft, bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
	long charsnow;
#endif
	char **s, **p;
	char *t;

	s = (char **)evalstrs;
	while( n-- > 0){
		if( *s == 0) break;
		while( *s++ != 0){}
	}
	if( *s == 0) return 0;
	bytesleft = 0;
	charsleft = 0;
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft += bytesnow;
		charsleft += bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
		t = *p;
		charsleft -= WSCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
#endif
		++p;
	}


	WSPutNext( wslp, WSTKSTR);
#if EVALSTRS_AS_BYTESTRINGS
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		WSPut8BitCharacters( wslp, bytesleft, (unsigned char*)*p, bytesnow);
		++p;
	}
#else
	WSPut7BitCount( wslp, charsleft, bytesleft);
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		t = *p;
		charsnow = bytesnow - WSCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
		charsleft -= charsnow;
		WSPut7BitCharacters(  wslp, charsleft, *p, bytesnow, charsnow);
		++p;
	}
#endif
	return WSError( wslp) == WSEOK;
}
#endif /* CARDOF_EVALSTRS */


static int  _definepattern( WSLINK wslp, char* patt, char* args, int func_n)
{
	WSPutFunction( wslp, "DefineExternal", (long)3);
	  WSPutString( wslp, patt);
	  WSPutString( wslp, args);
	  WSPutInteger( wslp, func_n);
	return !WSError(wslp);
} /* _definepattern */


int _PREPDoCallPacket( WSLINK wslp, struct func functable[], int nfuncs)
{
	int len;
	int n, res = 0;
	struct func* funcp;

	if( ! WSGetInteger( wslp, &n) ||  n < 0 ||  n >= nfuncs) goto L0;
	funcp = &functable[n];

	if( funcp->f_nargs >= 0
	&& ( ! WSTestHead(wslp, "List", &len)
	     || ( !funcp->manual && (len != funcp->f_nargs))
	     || (  funcp->manual && (len <  funcp->f_nargs))
	   )
	) goto L0;

	stdlink = wslp;
	res = (*funcp->f_func)( wslp);

L0:	if( res == 0)
		res = WSClearError( wslp) && WSPutSymbol( wslp, "$Failed");
	return res && WSEndPacket( wslp) && WSNewPacket( wslp);
} /* _PREPDoCallPacket */


wsapi_packet PREPAnswer( WSLINK wslp)
{
	wsapi_packet pkt = 0;
	int waitResult;

	while( ! PREPDone && ! WSError(wslp)
		&& (waitResult = WSWaitForLinkActivity(wslp),waitResult) &&
		waitResult == WSWAITSUCCESS && (pkt = WSNextPacket(wslp), pkt) &&
		pkt == CALLPKT)
	{
		PREPAbort = 0;
		if(! PREPDoCallPacket(wslp))
			pkt = 0;
	}
	PREPAbort = 0;
	return pkt;
} /* PREPAnswer */



/*
	Module[ { me = $ParentLink},
		$ParentLink = contents of RESUMEPKT;
		Message[ MessageName[$ParentLink, "notfe"], me];
		me]
*/

static int refuse_to_be_a_frontend( WSLINK wslp)
{
	int pkt;

	WSPutFunction( wslp, "EvaluatePacket", 1);
	  WSPutFunction( wslp, "Module", 2);
	    WSPutFunction( wslp, "List", 1);
		  WSPutFunction( wslp, "Set", 2);
		    WSPutSymbol( wslp, "me");
	        WSPutSymbol( wslp, "$ParentLink");
	  WSPutFunction( wslp, "CompoundExpression", 3);
	    WSPutFunction( wslp, "Set", 2);
	      WSPutSymbol( wslp, "$ParentLink");
	      WSTransferExpression( wslp, wslp);
	    WSPutFunction( wslp, "Message", 2);
	      WSPutFunction( wslp, "MessageName", 2);
	        WSPutSymbol( wslp, "$ParentLink");
	        WSPutString( wslp, "notfe");
	      WSPutSymbol( wslp, "me");
	    WSPutSymbol( wslp, "me");
	WSEndPacket( wslp);

	while( (pkt = WSNextPacket( wslp), pkt) && pkt != SUSPENDPKT)
		WSNewPacket( wslp);
	WSNewPacket( wslp);
	return WSError( wslp) == WSEOK;
}


int PREPEvaluate( WSLINK wslp, char* s)
{
	if( PREPAbort) return 0;
	return WSPutFunction( wslp, "EvaluatePacket", 1L)
		&& WSPutFunction( wslp, "ToExpression", 1L)
		&& WSPutString( wslp, s)
		&& WSEndPacket( wslp);
} /* PREPEvaluate */


int PREPEvaluateString( WSLINK wslp, char* s)
{
	int pkt;
	if( PREPAbort) return 0;
	if( PREPEvaluate( wslp, s)){
		while( (pkt = PREPAnswer( wslp), pkt) && pkt != RETURNPKT)
			WSNewPacket( wslp);
		WSNewPacket( wslp);
	}
	return WSError( wslp) == WSEOK;
} /* PREPEvaluateString */


void PREPDefaultHandler( WSLINK wslp, int message, int n)
{
	wslp = (WSLINK)0; /* suppress unused warning */
	n = 0; /* suppress unused warning */

	switch (message){
	case WSTerminateMessage:
		PREPDone = 1;
	case WSInterruptMessage:
	case WSAbortMessage:
		PREPAbort = 1;
	default:
		return;
	}
}


static int _WSMain( char **argv, char **argv_end, char *commandline)
{
	WSLINK wslp;
	int err;

	if( !stdenv)
		stdenv = WSInitialize( (WSEnvironmentParameter)0);
	if( stdenv == (WSEnvironment)0) goto R0;

	if( !stdhandler)
		stdhandler = (WSMessageHandlerObject)PREPDefaultHandler;
	wslp = commandline
		? WSOpenString( stdenv, commandline, &err)
		: WSOpenArgv( stdenv, argv, argv_end, &err);

	if( wslp == (WSLINK)0){
		WSAlert( stdenv, WSErrorString( stdenv, err));
		goto R1;
	}

	if( stdyielder) WSSetYieldFunction( wslp, stdyielder);
	if( stdhandler) WSSetMessageHandler( wslp, stdhandler);

	if( WSInstall( wslp))
		while( PREPAnswer( wslp) == RESUMEPKT){
			if( ! refuse_to_be_a_frontend( wslp)) break;
		}

	WSClose( wslp);
R1:	WSDeinitialize( stdenv);
	stdenv = (WSEnvironment)0;
R0:	return !PREPDone;
} /* _WSMain */


int WSMainString( char *commandline)
{
	return _WSMain( (char **)0, (char **)0, commandline);
}

int WSMainArgv( char** argv, char** argv_end) /* note not FAR pointers */
{   
	static char FAR * far_argv[128];
	int count = 0;
	
	while(argv < argv_end)
		far_argv[count++] = *argv++;
		 
	return _WSMain( far_argv, far_argv + count, (charp_ct)0);

}

int WSMain( int argc, char ** argv)
{
 	return _WSMain( argv, argv + argc, (char *)0);
}
