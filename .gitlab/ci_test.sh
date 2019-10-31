$CI_PROJECT_DIR/build/test_rutherford_inspect \
-n 1000 \
-o /data/test_rutherford_inspect.root \
-d 1 \
-v 10 \
-S 10

python3 $CI_PROJECT_DIR/.gitlab/upload_qa \
https://gitlab.cern.ch \
mkwon/RooParticle \
$GITLAB_API_KEY \
/data/test_rutherford_inspect.root
$CI_COMMIT_SHA
