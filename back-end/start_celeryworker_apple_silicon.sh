# `--pool solo` is a workaround to ensure that this works on apple silicon. Since we only use one worker anyway, this is
#  fine (https://github.com/celery/celery/issues/7324#issuecomment-1990191006)
celery -A tasks worker --loglevel=INFO --pool solo