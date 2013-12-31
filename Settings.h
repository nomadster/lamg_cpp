/*
 * File:   Settings.h
 * Author: nomadsoul
 *
 * Created on 5 settembre 2013, 11.21
 */

#ifndef SETTINGS_H
#define	SETTINGS_H

/** In C++ è brutto usare le macro, ci sarebbero gli enum ei const
 *  Per adesso me ne frego.
 * TODO: Cambiare le macro con enum e/o const.
 * TODO2: Rendere questi poi parametri in un file dei parametri, che vengono
 * parsati all'avvio del programma.
 */

/**************************
 * SETUP - Par. Generali  *
 **************************/
/** Determina se il multigrid viene costruito in SetD() utilizzando
 *  il peso degli archi [FALSE] oppure in SetGraph() utilizzando la sola
 *  struttura del grafo [TRUE].
 */
#define SETUP_BUILD_STATIC_MULTIGRID true
/*
minCoarseSize = 300 %200                % A guide for the minimum coarsest level size
        maxCoarseRelaxAcf = 0.3 %0.7            % If relaxation converges at this rate or faster, this becomes the coarsest level
        relaxAcfMinSweeps = 7                   % Minimum number of sweeps to run to estimate relax ACF
        setupNumAggLevels = 100                 % maximum #aggregation coarsening levels to construct during the setup phase
        setupNumLevels = 100                    % maximum TOTAL # of coarsening levels to construct during the setup phase
        setupSave = false                       % If true, saves level hierarchy to an m-file                %
        tvNumLocalSweeps = 0                    % # local sweeps to perform on each TV at each node i during coarsening (then retract)
        tvInitialGuess = 'random'               % Type of TV initial guess (random/geometric)
        interpType = 'caliber1'                 % Interpolation operator type
        restrictType = 'transpose'              % Restriction operator type
        nuDesign = 'split_evenly_post'          % Strategy of splitting the relaxation sweep number nu at each level to nuPre and nuPost
 */

// #initial test vectors (TVs) at each level.
#define TV_NUM 4
// #TVs to add upon each aggregation coarsening
#define TV_INC 1
// Maximum allowed#TVs
#define TV_MAX 10
// #global sweeps to perform on each initial TV
#define SETUP_TV_SWEEPS 4
// Max size for direct solver
#define MAX_DIRECT_SOLVE_SIZE 200
// Maximum #aggregation coarsening levels to construct during the setup phase
#define SETUP_MAX_AGG_LEVELS 100
// Maximum TOTAL # of coarsening levels to construct during the setup phase
#define SETUP_MAX_LEVELS 100
// Solution cycle index. Also the design cycle index during setup phase.
#define SETUP_CYCLE_INDEX 1.5
// Minimum number of sweeps to run to estimate relax ACF
#define SETUP_RELAX_ACF_MIN_SWEEPS 7
//If relaxation converges at this rate or faster, this becomes the CoarsestLevel
#define SETUP_MAX_COARSE_RELAX_ACF 0.3

/**************************
 * SETUP - Elimination    *
 **************************/

// Se falso vengono usati solo livelli di tipo AGGREGATION
//UNUSED->#define SETUP_ELIMINATION true
// Il grado massimo dei nodi low-degree
#define SETUP_ELIMINATION_MAX_DEGREE 4
//Numero massimo di stage di eliminazione
#define SETUP_ELIMINATION_MAX_STAGES 1000
/* Frazione minima di nodi eliminati richiesta nella fase di ELIMINATION
 * prima che si passi ad una fase di AGGREGATION */
#define SETUP_ELIMINATION_MIN_ELIM_FRACTION 0.01

/**************************
 * SETUP - Aggregation    *
 **************************/

#define SETUP_AGGREGATION_WEAK_EDJE_THRESHOLD 0.1
// % Mark all locally-high-degree nodes as seeds
//t = obj.options.aggregationDegreeThreshold;
#define SETUP_AGGREGATION_DEGREE_THRESHOLD 8
// % #sweeps (nu) to use during coarsening stages.
#define SETUP_NU_DEAFULT 3
//% Bound on gamma*alpha during coarsening ???
#define SETUP_COARSENING_WORK_GUARD 0.7

// Il nostro alpha_max
// NON USO PIU' #define MAX_COARSENING_RATIO ((SETUP_COARSENING_WORK_GUARD)/(SETUP_CYCLE_INDEX))

//min # coarsening stages to generate
#define SETUP_MIN_AGGREGATION_STAGES 1
//max # coarsening stages to generate before deciding on the best
#define SETUP_MAX_AGGREGATION_STAGES 2
//Quante passate di GS fare al massimo
#define SOLVE_MAX_COARSEST_SWEEP 400
//Ogni quanto controllo la convergenza
#define SOLVE_RELAX_FREQUENCY 10

// Maximum energy ratio in delta model
#define RATIO_MAX 2.5

#endif	/* SETTINGS_H */
