import logging

from block_execution import execute_blocks
from celery import Celery

# We keep celery in its own file, so we can run without needing celery or redis whilst running tests.

logger = logging.getLogger()

logger.info('Intialising celery')
celery = Celery(__name__, broker="redis://localhost:6379/0")
celery.conf.result_backend = "redis://localhost:6379/0"
celery.conf.broker_transport_options = {'visibility_timeout': 20}  # 20 seconds
celery.conf.task_serializer = 'json'


@celery.task()
def execute_blocks_celery(message, session_id, socketio_url):
    execute_blocks(message, session_id, socketio_url)
