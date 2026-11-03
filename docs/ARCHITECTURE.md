# EES Architecture Document
## Engineering Equation Solver — Numerical Methods & Engineering Analysis Workbench

---

## 1. System Overview & Philosophy

**EES (Engineering Equation Solver)** is an integrated scientific computing workbench and educational laboratory engineered for Computational Methods / Numerical Methods curriculum and engineering problem solving.

### Architectural Philosophy
```text
Solve → Explain → Visualize → Compare → Understand
```
1. **Zero External Numerical Dependencies**: All algorithms (roots, optimization, matrices, eigenvalues, splines, ODEs, PDEs) are implemented from first principles in C++ (C++17/20).
2. **Layered Separation of Concerns**: Mathematical foundations (`Vector`, `Matrix`, `ExpressionParser`, `ErrorAnalysis`, `NumericalResult`) are strictly isolated from control flow and user interface logic.
3. **Structured Introspection**: Every solver records complete iteration trajectories ($x_k, f(x_k), \epsilon_k$, auxiliary variables, convergence metrics) to support interactive step-by-step examination, error graphing, and automated multi-method comparison.
4. **Cinematic & Restrained Engineering UX**: Clean visual hierarchy, ANSI/terminal scientific display and web/GUI rendering capabilities, monochrome dark graphite aesthetic with restrained color accents (amber for active computation, green for convergence, red for failure).

---

## 2. Layered Architecture

```text
┌────────────────────────────────────────────────────────────────────────┐
│                        User Interface Layer                            │
│  - Interactive CLI Workbench (Menu, Step-by-Step, Compare, Benchmarks) │
│  - Rich Text / Terminal ASCII & Unicode Canvas Visualizations          │
│  - HTML5/Canvas Exportable Visualizer Dashboard                       │
└────────────────────────────────────┬───────────────────────────────────┘
                                     │
┌────────────────────────────────────▼───────────────────────────────────┐
│                      Application & Controller Layer                    │
│  - Problem Library Manager (Mechanical, Civil, Electrical, Thermal)    │
│  - Numerical Method Advisor (Heuristic Recommendation Engine)          │
│  - Method Comparison Engine (Side-by-side performance & error matrix)  │
│  - Interactive Stepper / Debugger Controller                           │
└────────────────────────────────────┬───────────────────────────────────┘
                                     │
┌────────────────────────────────────▼───────────────────────────────────┐
│                   Numerical Methods Engine (Syllabus)                  │
│  ┌───────────────────────┬───────────────────────┐                     │
│  │ UNIT I: Foundations,  │ UNIT II:              │                     │
│  │ Roots & Optimization  │ Interpolation & Quad  │                     │
│  ├───────────────────────┼───────────────────────┤                     │
│  │ UNIT III:             │ UNIT IV:              │                     │
│  │ Linear Systems &      │ ODE Solvers &         │                     │
│  │ Matrix Factorization  │ PDE Solvers           │                     │
│  └───────────────────────┴───────────────────────┘                     │
└────────────────────────────────────┬───────────────────────────────────┘
                                     │
┌────────────────────────────────────▼───────────────────────────────────┐
│                   Core Mathematical Infrastructure                     │
│  - Vector (Dynamic, L1/L2/Linf norms, dot/cross, linspace)            │
│  - Matrix (Row-major, LU det, Elementary Row Ops, norms, transpose)   │
│  - ExpressionParser (AST / Shunting-yard RPN, safe evaluation, diff)   │
│  - ErrorAnalysis (Absolute, Relative, %, IEEE-754, Catastrophic cancel)│
│  - NumericalResult & IterationRecord (Structured trajectory models)   │
└────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Directory & Module Organization

```text
EES/
├── CMakeLists.txt              # Unified build configuration (Ninja, GCC, MSVC, Clang)
├── README.md                   # Top-level showcase and syllabus mapping
├── docs/
│   ├── ARCHITECTURE.md         # Detailed system architecture (this document)
│   ├── ALGORITHMS.md           # Formulae, pseudocode, stability, complexity
│   ├── MATHEMATICAL_FOUNDATION.md # Numerical analysis theory
│   ├── TESTING.md              # Test strategy and benchmark problems
│   ├── USER_GUIDE.md           # Operational manual and tutorials
│   └── SYLLABUS_MAPPING.md     # Direct syllabus topic-to-file mapping
├── include/ees/
│   ├── core/                   # Vector, Matrix, ExpressionParser, ErrorAnalysis, Types
│   ├── unit1/                  # Taylor, Foundations, Roots (Bisection, Newton, Secant), Opt
│   ├── unit2/                  # Interpolation (Diff tables, Newton, Lagrange), Integration
│   ├── unit3/                  # Linear Systems (Gauss, Jordan, LU, Cholesky), Splines, Power
│   ├── unit4/                  # ODEs (Euler, RK4, Predictor-Corrector) & PDEs (Heat, Wave, Laplace)
│   ├── comparison/             # Multi-method benchmark and comparison engine
│   ├── advisor/                # Method recommendation heuristic system
│   ├── problems/               # Domain-specific engineering problem catalog
│   └── ui/                     # Terminal graphics, charts, tables, interactive shell
├── src/                        # Implementations matching include/ees/
├── tests/                      # Comprehensive automated test suites per unit
├── examples/                   # Standalone reproducible tutorial scripts
└── build/                      # Build output directory
```

---

## 4. Key Design Decisions & Resolved Ambiguities

1. **Memory Representation**:
   - `Matrix` uses a contiguous 1D row-major array (`std::vector<double>`) for maximum cache locality and SIMD vectorization potential.
   - Matrix row operations (`swapRows`, `addRowMultiple`) modify the matrix directly, avoiding excessive copies during Gaussian elimination.
2. **Expression Evaluation & Differentiation**:
   - Built with an internal token-based Shunting-Yard RPN parser to guarantee zero external dependencies.
   - Supports analytic-grade numerical differentiation using 2nd and 4th order central finite differences ($O(h^2)$ and $O(h^4)$ error), enabling Newton-Raphson, Secant, and Optimization methods to operate seamlessly without manual derivative entry.
3. **Step-by-Step State Machines**:
   - Every iterative method is written with an internal generator or step-aware loop that can either run to completion or advance one iteration at a time for interactive student viva demonstrations.
4. **Visualization Strategy**:
   - **Terminal ANSI/Unicode Canvas**: Generates high-density terminal scatter plots, contour slices, interval shrinkage views, and difference tables that render directly in Windows PowerShell / Command Prompt.
   - **Data Export**: Generates clean JSON and tabular data suitable for graphing or downstream analysis.

---

## 5. Development Phases Roadmap

- [x] **Phase 1: Project Architecture & Core Infrastructure** (Vector, Matrix, Parser, ErrorAnalysis, CMake, Ninja, Core Tests)
- [ ] **Phase 2: Unit I — Numerical Foundations, Computer Arithmetic & Taylor Series**
- [ ] **Phase 3: Unit I — Root-Finding Laboratory & Method Comparison**
- [ ] **Phase 4: Unit I — Unconstrained & Multivariate Optimization**
- [ ] **Phase 5: Unit II — Interpolation & Finite Differences Laboratory**
- [ ] **Phase 6: Unit II — Numerical Integration & Quadrature Laboratory**
- [ ] **Phase 7: Unit III — Linear Systems, Pivoting, LU/Cholesky & Power Method**
- [ ] **Phase 8: Unit III — Spline Approximation Laboratory**
- [ ] **Phase 9: Unit IV — Numerical Solution of ODEs**
- [ ] **Phase 10: Unit IV — Numerical Solution of PDEs**
- [ ] **Phase 11: Visualization Engine & Terminal Graphics Canvas**
- [ ] **Phase 12: Method Comparison Engine & Engineering Problem Library**
- [ ] **Phase 13: Numerical Method Advisor**
- [ ] **Phase 14: Comprehensive Automated Test Suite**
- [ ] **Phase 15: Full Documentation Suite (Syllabus Mapping, Algorithms, Guide)**
- [ ] **Phase 16: Interactive UI Polish & PBL Viva Demonstration Mode**

