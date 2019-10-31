echo "Rutherford Scattering Simulation"

$CI_PROJECT_DIR/build/test_rutherford_inspect \
-n 50 \
-o /data/test_rutherford_inspect.root \
-d 1 \
-v 10 \
-S 10

echo "Upload to Gitlab Repository Board"

python3 $CI_PROJECT_DIR/.gitlab/upload_qa.py \
https://gitlab.cern.ch \
mkwon/RooParticle \
$GITLAB_API_KEY \
/data/test_rutherford_inspect.root \
$CI_COMMIT_SHA
