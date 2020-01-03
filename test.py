
import os,logging

# logging
logger=logging.getLogger('AIM-UP')
logger.setLevel(logging.DEBUG)

ch=logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter=logging.Formatter('%(asctime)s/%(name)s/%(levelname)s/%(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

ch=logging.FileHandler('log.txt')
formatter=logging.Formatter('%(asctime)s/%(name)s/%(levelname)s/%(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

logger.debug('debug')
logger.info('info')
logger.warning('This is a warning')
logger.error('This is an error')
logger.critical('This is an critical')
