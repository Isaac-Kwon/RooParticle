# 남은 작업

2026-11-13 기준. 사용법은 [README.md](README.md)와 Doxygen 문서(`make doc`)에 있다.
이 프로젝트는 3체 이상의 산란을 전제로 한다. 2체 닫힌 해는 검증 기준과 특수 경우로만 쓴다.

## 2026-11-13에 한 것

빌드 35개 타깃, `ctest` 21개 통과, `make doc` 경고 0건 (WSL, ROOT 6.36.04).

- **적분기**: `TAYLOR`(기본), `SEULER`, `VERLET`, `RK4`, 가변 스텝 `RK45` 중 선택. 상대론 모드(운동량 적분).
  작용·반작용으로 힘 계산 절반화. 운동량·각운동량 조회.
- **종료 조건**: OR 결합, abort inspector, 스텝 상한, 종료 사유 조회, 메서드 enum화, `atan2` 각도.
  새 메서드 16종. N체용 `ASYD`/`ASYV`(남은 편향량 기준) 포함.
- **점근 상태**: 2체 정확해 `asymptote`, N체용 `asymptoteN`(꼬리 보정, 쌍별 합, 수치 기준값, 묶인 계 감지, 오차 추정).
- **힘**: 중력 구현, Yukawa, Woods-Saxon, softening. 쿨롱 상수·중력 상수 설정 가능.
- **궤적**: DCA 포물선 보간, 기록 간격.
- **테스트 프로그램**: 결과 저장 복구, 시드·작업 번호 통일, 실행 메타데이터 기록, 스텝 상한.
- **빌드**: 정적 라이브러리, 테스트 자동 등록, GitHub Actions 워크플로, `.gitignore`.
- **스크립트**: `runmulti.py`·`runmulti.sh` 재작성, `plot_rutherford.C`, `even_odd` 매크로 인자화.
- **문서**: Doxygen 설정, 전 헤더 API 주석, 튜토리얼·개념 설명·검증 결과 페이지.
- **검증**: 러더퍼드 산란각, 최근접거리, 표적 반동, 미분 단면적, 상대론, 핵력, 3·4체 보존량과 수렴.

고친 버그:

- 고정점 DCA가 항상 0이던 문제 (`path::GetDCA`).
- `oddevent`에서 `Invincible` 입자가 스텝을 세지 않아 `CNT` 조건이 끝나지 않던 문제 (`main`이 끝나지 않던 원인).
- 복사된 `inspectorE`가 원본 이벤트를 감시하던 문제, 복사 시 `inspectorP`가 엉뚱한 입자에 연결될 수 있던 문제.
- `eventVGeneral(x, y, z)` 생성자가 박스 크기를 설정하지 않던 문제.
- `recorder::MinSave`가 저장하지 않던 문제 (최소 DCA 대신 초기값이 기록되고 있었다).
- 중력 상수의 지수 (1.2331e-17 → 1.2331e-39).
- inspector·`pEvent` 메모리 누수, 입자 복사 시 스텝 수 누락, 초기화 안 된 포인터들.
- RK45가 인력 특이점에서 한 호출 안에 끝나지 않던 문제 (종료 사유 `FAILED`).
- `main`·`test_collider_HeHe`와 `runmulti.py`의 인자 순서 불일치.
- `runmulti.py`가 끝나지 않던 문제, `test_collider_HeHe`의 인자 없는 실행.

## 설계 메모

- inspector의 AND 결합은 의도된 것이다. 입자를 최대한 멀리 보내기 위해
  "점근 각도"와 "최소 스텝 수"를 둘 다 만족해야 멈춘다. `CNT`는 상한이 아니라 하한이다.
- 다체에서는 출사 방향 오차가 정지 거리 R에 대해 1/R로만 줄어든다 (고정 표적 2체는 1/R²).
  1e-6 rad에 3체 자유는 약 1.3e6 fm가 필요하다. `asymptoteN`의 보정은 이 오차를 1/R²로 줄인다.

## 1. 결정이 필요한 것

- [ ] **기본 적분기를 `VERLET`으로 바꿀지**. 같은 비용에 훨씬 정확하지만 기존 결과와 수치가 달라진다.
- [ ] **쿨롱 상수**: 코드의 0.00154629는 1.4404 MeV·fm에 해당한다. CODATA 기준(1.439965 MeV·fm)으로는 0.00154587이고 0.03% 차이다.
  바꾸면 모든 결과 수치가 달라진다. 지금은 `SetK()`로만 바꿀 수 있다.
- [ ] **`Invincible` 입자의 시간**: 진행하지 않는다. 그래서 움직이는 입자의 궤적과 비교하는 `GetDCA(path*)`가
  시간 불일치로 -1을 돌려준다 (고정점 오버로드는 동작). 시간을 진행시킬지 결정 필요.
- [ ] **`test/main_rutherford_track.cpp`**: 컴파일되지 않아 빌드에서 제외돼 있다. 선언은 단일 충돌계수용,
  본문은 몬테카를로용이라 어느 쪽인지 정해야 고칠 수 있다.
- [ ] **LICENSE 종류**.
- [ ] **`SetInitial_D`, `preEvent`, `event::preDerive`·`SetInitial`**(초기 위치 자동 설정): 완성 또는 삭제.
  `event::SetInitial()`은 선언만 있어 호출하면 링크 에러다.
- [ ] 병합이 끝난 원격 브랜치 정리.

## 2. 확인이 필요한 것 (실행해 보지 못함)

- [ ] **GitHub Actions** (`.github/workflows/ci.yml`): YAML 문법만 확인했다. 푸시해서 `rootproject/root` 이미지에서 도는지 확인.
- [ ] **Condor 파일** (`macro/condor/`): 실행 파일 이름 오타를 고치고 스크립트에 `set -e`, `make -j4`를 넣었다. 클러스터에서 확인 필요.
  clone 주소가 아직 CERN GitLab이다. 또 `title ClusterId ProcId`를 넘기므로 작업 번호가 ClusterId, 시드가 ProcId가 되어
  `queue N`(N>1)이면 한 클러스터의 모든 작업이 같은 파일에 쓴다. `ProcId`를 작업 번호로 넘기도록 바꿔야 한다.
- [ ] `Inspect()`를 직접 도는 테스트 5개의 스텝 상한 경로 (상한이 1억 스텝이라 실제로 걸리는 경우를 돌려 보지 못했다).
- [ ] README 배지를 GitHub Actions(`Isaac-Kwon/RooParticle`)로 바꿨다. 워크플로가 한 번 돌기 전에는 상태가 표시되지 않는다.

## 3. 알려진 한계

- [ ] **Verlet 캐시**: 스텝 사이에 `SetK` 등 힘 파라미터나 속도·질량을 바꾸면 첫 스텝에 반영되지 않는다.
- [ ] **이벤트 복사**: 입자를 일반 `particle`로 복사한다 (전하는 유지). `new event(*oddevent)`는 일반 `event`가 된다.
  `experiment::MakeEvent`가 이 경로라 `oddevent` 템플릿을 쓸 수 없다.
- [ ] **`experiment`**: `MakeEvent`를 다시 부르면 이전 이벤트가 누수된다. 가변 인자의 개수·타입 검사가 없다.
- [ ] **`asymptote`(2체)**: 멀어지는 방사형 인력 쌍(E>0)에 점근 상태가 있는데도 없다고 답한다. `asymptoteN`은 자체 처리한다.
- [ ] **`asymptoteN`**: 오차 추정은 2~4체, 쿨롱, 쿨롱+Woods-Saxon, 상대론 1차에서 검증했다.
  중력·Yukawa·softening·5체 이상·`oddevent`는 미검증. 상대론 모드는 1차 보정만 된다.
  `coulombForce`를 상속한 클래스는 닫힌 식 대신 수치 적분 경로를 타서 호출당 약 60배 느리다 (결과는 같다).
  `ASYD`/`ASYV`는 각 입자 자신의 속도 기준이라, 느리게 반동하는 표적이 종료 시점을 결정한다.
- [ ] **RK45 실패 후 상태**: `FAILED`로 끝난 이벤트의 위치·속도는 NaN일 수 있다. 그 이벤트는 버려야 한다.
- [ ] **RK45**: `DeriveDT(dt)`가 내부 스텝을 `dt` 이하로 제한하므로 작은 `dt`로 부르면 허용 오차가 무의미하다.
  상대 허용 오차가 좌표 원점 기준이다.
- [ ] **상대론 모드**: 지연 효과와 자기력은 없다. `TAYLOR`와 `oddevent`는 설정을 무시한다.
- [ ] **핵력**: 이벤트에 넣으면 모든 입자쌍에 작용한다. 특정 쌍에만 거는 기능이 없다.
  Yukawa 기본 세기는 예시값이고, Woods-Saxon은 실험 데이터와 비교하지 않았다.
- [ ] `inspectorE`의 `DXL`·`DNG`는 `Evaluate()`가 돌려주는 값과 `Inspect()`가 비교하는 값이 다르다.
- [ ] `MAXSTEPS`는 상한이 아니라 `DeriveMAX`가 항상 도는 고정 스텝 수다 (이름이 오해를 부른다).
- [ ] `plot_rutherford.C`는 명목 충돌계수와 속도로 비교한다. 잔차(최대 0.024 rad)는 유한 시작 거리의 영향이다.
- [ ] `git_sha` 메타데이터는 커밋하지 않은 변경을 표시하지 않는다.

## 4. 더 할 수 있는 것

- [ ] 검증 추가: `oddevent`, Yukawa, 중력, 5체 이상.
- [ ] `test_integrator_compare`처럼 `ASYD`/`ASYV`·`FAILED`·상대론 inspector의 전용 테스트.
- [ ] 테스트마다 복붙된 `FixedTargetExperiment`/`ColliderExperiment`를 `experiment`로 통합.
- [ ] 분석 스크립트를 점근값(`asymptote`/`asymptoteN`) 기준으로.
- [ ] ROOT 딕셔너리 또는 PyROOT 바인딩.
- [ ] Graphviz 설치 시 Doxygen 클래스 다이어그램 (지금은 `dot`이 없어 꺼져 있다).
- [ ] 죽은 주석 코드 정리 (`eventViewer`, `preEvent`는 전체가 주석이다).
