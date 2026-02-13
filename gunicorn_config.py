bind = "127.0.0.1:8000"
workers = 1
worker_class = "sync"
timeout = 120
graceful_timeout = 30
accesslog = "/var/log/nanosim/access.log"
errorlog = "/var/log/nanosim/error.log"
loglevel = "info"

# Memory optimizations
max_requests = 100
max_requests_jitter = 20
worker_tmp_dir = "/dev/shm"
preload_app = True

# Connection handling
keepalive = 5
