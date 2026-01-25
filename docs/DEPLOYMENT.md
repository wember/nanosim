# Deployment Guide: AWS Lightsail

This guide walks through deploying the simulation plot browser to AWS Lightsail as a subdomain of your existing domain.

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Step 1: Create Lightsail Instance](#step-1-create-lightsail-instance)
- [Step 2: Configure DNS](#step-2-configure-dns)
- [Step 3: Configure Lightsail Firewall](#step-3-configure-lightsail-firewall)
- [Step 4: Connect and Setup Server](#step-4-connect-and-setup-server)
- [Step 5: Deploy Application](#step-5-deploy-application)
- [Step 6: Configure Gunicorn](#step-6-configure-gunicorn)
  - [Step 6b: Configure Swap Space](#step-6b-configure-swap-space-essential-for-512mb-instances)
  - [Step 6c: Network and System Optimizations](#step-6c-network-and-system-optimizations)
- [Step 7: Configure Nginx](#step-7-configure-nginx)
- [Step 8: Setup SSL with Let's Encrypt](#step-8-setup-ssl-with-lets-encrypt)
- [Step 9: Additional System Hardening](#step-9-additional-system-hardening)
- [Step 10: Setup Automated Health Monitoring](#step-10-setup-automated-health-monitoring)
- [Step 11: Configure Log Rotation](#step-11-configure-log-rotation)
- [Step 12: Upload Simulation Data](#step-12-upload-simulation-data)
- [Step 13: Verify Deployment](#step-13-verify-deployment)
- [Maintenance Tasks](#maintenance-tasks)
- [Optional: Add Basic Authentication](#optional-add-basic-authentication)
- [Troubleshooting](#troubleshooting)
- [Cost Breakdown](#cost-breakdown)
- [Scaling Up](#scaling-up)
- [Security Recommendations](#security-recommendations)
- [Next Steps](#next-steps)

## Overview

**Final Setup:**

- URL: `https://plots.yourdomain.com`
- Service: AWS Lightsail ($3.50-$5/month)
- Server: Ubuntu 22.04 with Python/Flask
- Web Server: Nginx + Gunicorn
- SSL: Let's Encrypt (free, auto-renewal)
- Storage: 20-40GB SSD (included)
- Memory: 512MB minimum (with swap) or 1GB recommended

## Prerequisites

- AWS account with billing enabled
- Your domain managed in Route 53 (or ability to add DNS records)
- SSH client installed locally
- Git installed locally

## Step 1: Create Lightsail Instance

1. **Log into AWS Lightsail**
   - Go to https://lightsail.aws.amazon.com/
   - Click "Create instance"

2. **Configure Instance**
   - **Location**: Choose closest region to you
   - **Platform**: Linux/Unix
   - **Blueprint**: OS Only → Ubuntu 22.04 LTS
   - **Plan**:
     - $3.50/month (512 MB RAM) - Minimum, requires swap setup
     - $5/month (1 GB RAM) - Recommended for better performance
     - Higher plans for larger datasets
   - **Instance name**: `nanosim-plots`
   - Click "Create instance"

3. **Wait for Instance to Start** (~2 minutes)

4. **Create Static IP**
   - Go to "Networking" tab
   - Click "Create static IP"
   - Attach to `nanosim-plots` instance
   - Name it: `nanosim-plots-ip`
   - Note the IP address (e.g., 54.123.45.67)

## Step 2: Configure DNS

### Option A: Using Route 53 (Recommended)

1. **Go to Route 53**
   - Open AWS Console → Route 53
   - Select your hosted zone

2. **Create A Record**
   - Click "Create record"
   - Record name: `plots`
   - Record type: `A`
   - Value: Your Lightsail static IP (54.123.45.67)
   - TTL: 300 seconds
   - Click "Create records"

### Option B: Using External DNS Provider

Add an A record:

- Host: `plots` (or `plots.yourdomain.com`)
- Type: `A`
- Value: Your Lightsail static IP
- TTL: 300 or auto

**Wait 5-10 minutes for DNS propagation**, then verify:

```bash
nslookup plots.yourdomain.com
```

## Step 3: Configure Lightsail Firewall

1. **Go to Instance → Networking tab**
2. **Add firewall rules**:
   - SSH (TCP 22) - Already exists
   - HTTP (TCP 80) - Click "Add rule"
   - HTTPS (TCP 443) - Click "Add rule"

## Step 4: Connect and Setup Server

### SSH into Instance

**Via Lightsail Console:**

- Click "Connect using SSH" button

**OR via Terminal:**

```bash
# Download SSH key from Lightsail console first
chmod 400 ~/Downloads/LightsailDefaultKey-*.pem
ssh -i ~/Downloads/LightsailDefaultKey-*.pem ubuntu@plots.yourdomain.com
```

### Install Dependencies

```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Install Python and required packages
sudo apt install -y python3 python3-pip python3-venv nginx git

# Install supervisor for process management
sudo apt install -y supervisor

# Install certbot for SSL
sudo apt install -y certbot python3-certbot-nginx
```

## Step 5: Deploy Application

### Clone Repository

```bash
# Create app directory
sudo mkdir -p /var/www/nanosim
sudo chown ubuntu:ubuntu /var/www/nanosim

# Clone your repo (replace with your actual repo URL)
cd /var/www/nanosim
git clone https://github.com/wember/nanosim.git .

# OR copy files manually if not using git
# You can use SCP or Lightsail's file upload feature
```

### Setup Python Environment

```bash
cd /var/www/nanosim

# Create virtual environment
python3 -m venv venv

# Activate and install dependencies
source venv/bin/python
pip install --upgrade pip
pip install flask plotly numpy pandas scipy
```

### Create Data Directory

```bash
# Create data directory for simulation results
mkdir -p /var/www/nanosim/data

# Set permissions
chmod 755 /var/www/nanosim/data
```

## Step 6: Configure Gunicorn

### Create Gunicorn Configuration

```bash
sudo nano /var/www/nanosim/gunicorn_config.py
```

Add:

```python
bind = "127.0.0.1:8000"
workers = 1  # Use 1 worker for 512MB instance (2 for 1GB+)
worker_class = "sync"
timeout = 120
accesslog = "/var/log/nanosim/access.log"
errorlog = "/var/log/nanosim/error.log"
loglevel = "info"

# Memory optimizations
max_requests = 100  # Recycle workers to prevent memory leaks
max_requests_jitter = 20  # Randomize to avoid all workers restarting at once
worker_tmp_dir = "/dev/shm"  # Use RAM for temp files (faster)
preload_app = True  # Preload app before forking workers (saves memory)

# Connection handling
graceful_timeout = 30  # Time for graceful worker shutdown
keepalive = 5  # Keep connections alive for reuse
```

### Create Log Directory

```bash
sudo mkdir -p /var/log/nanosim
sudo chown ubuntu:ubuntu /var/log/nanosim
```

### Create Supervisor Configuration

```bash
sudo nano /etc/supervisor/conf.d/nanosim.conf
```

Add:

```ini
[program:nanosim]
directory=/var/www/nanosim
command=/var/www/nanosim/venv/bin/gunicorn -c gunicorn_config.py tools.browse_plots:app
user=ubuntu
autostart=true
autorestart=true
stopasgroup=true
killasgroup=true
stderr_logfile=/var/log/nanosim/error.log
stdout_logfile=/var/log/nanosim/access.log
startretries=10
startsecs=5
stopwaitsecs=10
redirect_stderr=true
environment=PYTHONUNBUFFERED="1"
```

### Start Gunicorn

```bash
sudo supervisorctl reread
sudo supervisorctl update
sudo supervisorctl start nanosim
sudo supervisorctl status nanosim  # Should show RUNNING

# Enable supervisor to start on boot
sudo systemctl enable supervisor
```

## Step 6b: Configure Swap Space (512MB Instance)

**Important for stability on small instances:**

```bash
# Create 1GB swap file
sudo fallocate -l 1G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

# Make swap permanent
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab

# Optimize swap usage (prefer RAM over swap)
sudo bash -c 'echo "vm.swappiness=10" >> /etc/sysctl.conf'
sudo sysctl vm.swappiness=10

# Verify swap is active
free -h
```

**Why this matters:**

- Prevents OOM (Out of Memory) kills
- Allows graceful degradation under memory pressure
- Essential for 512MB instances running plot generation
- Completely free (uses existing disk space)

## Step 6c: Network and System Optimizations

**Optimize system for better stability and performance:**

```bash
# Add network tuning
sudo tee -a /etc/sysctl.conf > /dev/null << 'EOF'

# Network optimizations for stability
net.core.somaxconn = 1024
net.ipv4.tcp_max_syn_backlog = 2048
net.ipv4.ip_local_port_range = 10000 65000
net.ipv4.tcp_tw_reuse = 1
net.ipv4.tcp_fin_timeout = 15

# TCP keepalive for detecting dead connections
net.ipv4.tcp_keepalive_time = 300
net.ipv4.tcp_keepalive_intvl = 30
net.ipv4.tcp_keepalive_probes = 5
net.ipv4.tcp_slow_start_after_idle = 0
net.core.netdev_max_backlog = 2500

# Memory/cache optimizations for low-memory systems
vm.vfs_cache_pressure = 50
vm.dirty_background_ratio = 5
vm.dirty_ratio = 10
EOF

# Apply changes
sudo sysctl -p
```

**What these do:**

- Increased connection queues for handling traffic bursts
- More ephemeral ports available
- Faster connection cleanup and reuse
- Better handling of concurrent connections
- TCP keepalive detects and closes dead connections (5 min timeout)
- Optimized memory cache behavior for limited RAM
- Better disk write performance under memory pressure

## Step 7: Configure Nginx

### Create Nginx Configuration

```bash
sudo nano /etc/nginx/sites-available/nanosim
```

Add:

```nginx
server {
    listen 80;
    server_name plots.yourdomain.com;

    client_max_body_size 100M;

    # Gzip compression for faster page loads
    gzip on;
    gzip_vary on;
    gzip_proxied any;
    gzip_comp_level 6;
    gzip_types text/plain text/css text/xml text/javascript application/json application/javascript application/xml+rss;

    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 300s;
        proxy_connect_timeout 75s;

        # Prevent request buffer overflow
        proxy_buffering on;
        proxy_buffer_size 4k;
        proxy_buffers 8 4k;
        proxy_busy_buffers_size 8k;

        # Retry on failure
        proxy_next_upstream error timeout http_502 http_503 http_504;
        proxy_next_upstream_tries 2;
    }

    # Security headers
    add_header X-Frame-Options "SAMEORIGIN" always;
    add_header X-Content-Type-Options "nosniff" always;
    add_header X-XSS-Protection "1; mode=block" always;

    # Optional: Add basic authentication
    # auth_basic "Restricted Access";
    # auth_basic_user_file /etc/nginx/.htpasswd;
}
```

### Enable Site

```bash
# Create symlink
sudo ln -s /etc/nginx/sites-available/nanosim /etc/nginx/sites-enabled/

# Remove default site
sudo rm /etc/nginx/sites-enabled/default

# Test configuration
sudo nginx -t

# Restart nginx
sudo systemctl restart nginx
```

### Test HTTP Access

Visit `http://plots.yourdomain.com` - you should see your app!

## Step 8: Setup SSL with Let's Encrypt

```bash
# Run certbot
sudo certbot --nginx -d plots.yourdomain.com

# Follow prompts:
# - Enter email address
# - Agree to terms
# - Choose to redirect HTTP to HTTPS (recommended)
```

Certbot will:

- Obtain SSL certificate
- Automatically update Nginx configuration to use HTTP/2
- Setup auto-renewal

**Important:** After certbot runs, verify HTTP/2 is enabled:

```bash
# Check nginx config
sudo grep "listen 443" /etc/nginx/sites-available/nanosim
# Should show: listen 443 ssl http2;
```

### Verify SSL Auto-Renewal

```bash
# Test renewal
sudo certbot renew --dry-run

# Check renewal timer
sudo systemctl status certbot.timer
```

## Step 9: Additional System Hardening

### Limit Systemd Journal Size

Prevent journal logs from consuming excessive disk space:

```bash
sudo mkdir -p /etc/systemd/journald.conf.d
sudo tee /etc/systemd/journald.conf.d/size-limit.conf > /dev/null << 'EOF'
[Journal]
SystemMaxUse=50M
SystemMaxFileSize=10M
MaxRetentionSec=7day
EOF

sudo systemctl restart systemd-journald
```

### Increase File Descriptor Limits

Prevent "too many open files" errors under load:

```bash
sudo tee -a /etc/security/limits.conf > /dev/null << 'EOF'

# Increase file descriptor limits for ubuntu user
ubuntu soft nofile 4096
ubuntu hard nofile 8192
EOF
```

### Create Daily Cleanup Script

Automatically clean temporary files and caches:

```bash
sudo tee /etc/cron.daily/nanosim-cleanup > /dev/null << 'EOF'
#!/bin/bash
# Clean up temporary files and optimize disk usage

# Remove old Python cache files
find /var/www/nanosim -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null

# Clean old pip cache
rm -rf /home/ubuntu/.cache/pip/* 2>/dev/null

# Vacuum systemd journal
journalctl --vacuum-time=7d >/dev/null 2>&1

# Clean apt cache
apt-get clean >/dev/null 2>&1
EOF

sudo chmod +x /etc/cron.daily/nanosim-cleanup
```

### Add Monitoring Aliases (Optional)

Add convenient commands for quick system checks:

```bash
cat >> ~/.bashrc << 'EOF'

# Aliases for monitoring
alias nanostatus="sudo supervisorctl status nanosim && free -h && df -h /"
alias nanologs="sudo tail -n 50 /var/log/nanosim/error.log"
alias nanohealth="ps aux | grep -E \"gunicorn|supervisor\" | grep -v grep; echo \"---\"; curl -s -o /dev/null -w \"HTTP %{http_code} - Time: %{time_total}s\n\" http://127.0.0.1:8000/"
EOF

source ~/.bashrc
```

## Step 10: Setup Automated Health Monitoring

Configure automated health checks to restart the application if it becomes unresponsive:

```bash
# Create health check cron job
sudo nano /etc/cron.d/nanosim-health
```

Add this configuration:

```cron
# Health check every 5 minutes
*/5 * * * * root curl -f http://localhost:5000/ > /dev/null 2>&1 || supervisorctl restart nanosim
```

Save and verify:

```bash
# Check cron job is registered
sudo ls -la /etc/cron.d/nanosim-health

# Verify syntax
sudo cat /etc/cron.d/nanosim-health

# Watch for health check activity (optional)
tail -f /var/log/syslog | grep CRON
```

**How it works:**

- Checks if application responds on port 8000 every 5 minutes
- If check fails (curl returns non-zero), supervisor restarts the app
- Provides self-healing capability for unexpected failures

## Step 11: Configure Log Rotation

Prevent log files from consuming disk space:

```bash
# Create logrotate configuration
sudo nano /etc/logrotate.d/nanosim
```

Add this configuration:

```logrotate
/var/log/nanosim/*.log {
    daily
    rotate 7
    compress
    delaycompress
    missingok
    notifempty
    create 0640 www-data www-data
    sharedscripts
    postrotate
        supervisorctl restart nanosim > /dev/null
    endscript
}
```

Test the configuration:

```bash
# Test logrotate (dry run)
sudo logrotate -d /etc/logrotate.d/nanosim

# Force rotation to verify it works
sudo logrotate -f /etc/logrotate.d/nanosim

# Verify rotated logs
ls -lh /var/log/nanosim/
```

**Configuration explained:**

- `daily`: Rotate logs every day
- `rotate 7`: Keep 7 days of logs
- `compress`: Compress old logs with gzip
- `delaycompress`: Don't compress most recent rotation (easier debugging)
- `missingok`: Don't error if log file is missing
- `notifempty`: Skip rotation if log is empty
- `postrotate`: Restart app after rotation to release file handles

## Step 12: Upload Simulation Data

### Option A: SCP from Local Machine

```bash
# From your local machine
scp -i ~/Downloads/LightsailDefaultKey-*.pem -r ./data/* ubuntu@plots.yourdomain.com:/var/www/nanosim/data/
```

### Option B: Use Lightsail File Upload

1. Go to Lightsail console
2. Click instance → Connect → Upload files
3. Upload data directory contents

### Option C: Git LFS (if data in repo)

```bash
# On server
cd /var/www/nanosim
git lfs pull
```

## Step 13: Verify Deployment

1. **Visit**: `https://plots.yourdomain.com`
2. **Check logs**:

   ```bash
   # Application logs
   tail -f /var/log/nanosim/access.log
   tail -f /var/log/nanosim/error.log

   # Nginx logs
   sudo tail -f /var/nginx/access.log
   sudo tail -f /var/nginx/error.log
   ```

3. **Verify all optimizations**:

   ```bash
   # Check swap is active
   free -h
   swapon --show

   # Verify network tuning (should show all optimizations)
   sudo sysctl net.core.somaxconn net.ipv4.tcp_max_syn_backlog net.ipv4.tcp_keepalive_time vm.vfs_cache_pressure

   # Check gunicorn worker count and config
   ps aux | grep gunicorn
   cat /var/www/nanosim/gunicorn_config.py | grep -E 'workers|graceful_timeout|keepalive'

   # Verify health check cron
   sudo cat /etc/cron.d/nanosim-health

   # Check log rotation config
   sudo cat /etc/logrotate.d/nanosim

   # Verify journal size limits
   sudo journalctl --disk-usage  # Should be under 50M

   # Check file descriptor limits
   ulimit -n  # Should be 4096

   # Test nginx gzip is enabled
   curl -H "Accept-Encoding: gzip" -I https://plots.yourdomain.com
   # Should see: Content-Encoding: gzip

   # Quick health check (if you added aliases)
   nanohealth
   ```

## Maintenance Tasks

### Monitor System Health

```bash
# Check memory usage (should show swap available)
free -h

# Monitor processes
htop

# Check application status
sudo supervisorctl status nanosim

# View recent health check activity
tail -20 /var/log/syslog | grep nanosim

# Check for OOM kills (should be none)
sudo dmesg | grep -i "killed process"
```

### Update Application Code

```bash
cd /var/www/nanosim
git pull
sudo supervisorctl restart nanosim
```

### Upload New Simulation Data

```bash
# From local machine
scp -i ~/key.pem -r ./data/20260122_* ubuntu@plots.yourdomain.com:/var/www/nanosim/data/
```

### View Logs

```bash
# Real-time logs
sudo supervisorctl tail -f nanosim

# All logs
sudo journalctl -u supervisor -f
```

### Restart Services

```bash
# Restart Flask app
sudo supervisorctl restart nanosim

# Restart Nginx
sudo systemctl restart nginx
```

### Check Disk Space

```bash
df -h
du -sh /var/www/nanosim/data/*
```

### Backup Data

```bash
# Create backup
cd /var/www/nanosim
tar -czf backup-$(date +%Y%m%d).tar.gz data/

# Download backup
scp -i ~/key.pem ubuntu@plots.yourdomain.com:/var/www/nanosim/backup-*.tar.gz ./
```

## Optional: Add Basic Authentication

To password-protect your site:

```bash
# Install apache2-utils
sudo apt install -y apache2-utils

# Create password file
sudo htpasswd -c /etc/nginx/.htpasswd yourusername

# Enter password when prompted

# Edit nginx config
sudo nano /etc/nginx/sites-available/nanosim

# Uncomment these lines:
# auth_basic "Restricted Access";
# auth_basic_user_file /etc/nginx/.htpasswd;

# Restart nginx
sudo systemctl restart nginx
```

## Troubleshooting

### App Won't Start

```bash
# Check supervisor status
sudo supervisorctl status nanosim

# View detailed logs
sudo supervisorctl tail nanosim stderr

# Restart supervisor
sudo systemctl restart supervisor
```

### 502 Bad Gateway

- Check if Gunicorn is running: `sudo supervisorctl status nanosim`
- Check Gunicorn logs: `tail -f /var/log/nanosim/error.log`
- Verify port 8000 is listening: `sudo netstat -tlnp | grep 8000`

### SSL Certificate Issues

```bash
# Check certificate status
sudo certbot certificates

# Force renewal
sudo certbot renew --force-renewal

# Check renewal logs
sudo cat /var/log/letsencrypt/letsencrypt.log
```

### Out of Disk Space

```bash
# Check space
df -h

# Find large directories
du -sh /var/www/nanosim/data/* | sort -h

# Clean old data
rm -rf /var/www/nanosim/data/20250101_*

# Clean logs
sudo journalctl --vacuum-time=7d
```

### Out of Memory / App Keeps Crashing

**Symptoms:**

- Service stops randomly
- `unix:///var/run/supervisor.sock no such file` errors
- OOM kills in logs: `sudo journalctl -u supervisor | grep OOM`

**Solutions:**

1. **Check if swap is configured:**

   ```bash
   free -h  # Should show swap space
   ```

2. **If no swap, add it (See Step 6b)**

3. **Reduce Gunicorn workers:**

   ```bash
   # Edit gunicorn config
   sudo nano /var/www/nanosim/gunicorn_config.py
   # Change workers = 2 to workers = 1

   # Restart
   sudo supervisorctl restart nanosim
   ```

4. **Check memory usage:**

   ```bash
   # Overall memory
   free -h

   # Top processes
   ps aux --sort=-%mem | head -10
   ```

5. **Monitor for crashes:**

   ```bash
   # Watch supervisor status
   watch -n 2 'sudo supervisorctl status nanosim'

   # Check for OOM kills
   sudo dmesg | grep -i kill

   # View recent health check activity
   tail -f /var/log/syslog | grep nanosim
   ```

6. **Verify all stability optimizations:**

   ```bash
   # Check swap
   swapon --show

   # Verify swappiness setting
   cat /proc/sys/vm/swappiness  # Should be 10

   # Check gunicorn config
   cat /var/www/nanosim/gunicorn_config.py | grep -E 'workers|max_requests|preload_app'

   # Verify network tuning
   sudo sysctl net.core.somaxconn net.ipv4.tcp_max_syn_backlog

   # Check health check is running
   sudo cat /etc/cron.d/nanosim-health
   ```

### Performance Issues / Slow Response

**Symptoms:**

- Site loads slowly
- Timeout errors
- High CPU usage

**Solutions:**

1. **Check nginx gzip compression:**

   ```bash
   # Test if gzip is working
   curl -H "Accept-Encoding: gzip" -I https://plots.yourdomain.com
   # Should see: Content-Encoding: gzip
   ```

2. **Monitor application performance:**

   ```bash
   # Check CPU/memory
   htop

   # Watch response times
   tail -f /var/log/nginx/access.log | grep --line-buffered -E 'GET|POST'
   ```

3. **Verify HTTP/2 is enabled:**

   ```bash
   # Check nginx SSL config
   sudo grep "listen 443" /etc/nginx/sites-available/nanosim
   # Should show: listen 443 ssl http2;
   ```

### Logs Growing Too Large

**Symptoms:**

- Disk space running low
- `/var/log/nanosim/` directory very large

**Solutions:**

1. **Verify log rotation is configured:**

   ```bash
   # Check config exists
   sudo cat /etc/logrotate.d/nanosim

   # Manually rotate logs
   sudo logrotate -f /etc/logrotate.d/nanosim

   # Check rotated logs
   ls -lh /var/log/nanosim/
   ```

2. **Clean old logs manually if needed:**

   ```bash
   # Remove old compressed logs
   sudo rm /var/log/nanosim/*.gz

   # Truncate current log files
   sudo truncate -s 0 /var/log/nanosim/*.log
   sudo supervisorctl restart nanosim
   ```

   ```

   ```

## Cost Breakdown

**512 MB Instance (Minimum):**

- **Lightsail Instance**: $3.50/month
- **Data Transfer**: Included (first 1 TB)
- **Static IP**: Included
- **SSL Certificate**: Free (Let's Encrypt)
- **DNS (Route 53)**: ~$0.50/month per hosted zone
- **Total: ~$4-5/month**

**1 GB Instance (Recommended):**

- **Lightsail Instance**: $5.00/month
- **Data Transfer**: Included (first 2 TB)
- **Static IP**: Included
- **SSL Certificate**: Free (Let's Encrypt)
- **DNS (Route 53)**: ~$0.50/month per hosted zone
- **Total: ~$5-6/month**

**Note:** 512MB instance requires swap configuration (Step 6b) for stability.

## Scaling Up

If you need more resources:

1. **Take Snapshot**
   - Lightsail console → Snapshots → Create snapshot

2. **Upgrade Instance**
   - Create new instance from snapshot
   - Choose larger plan ($10, $20, etc.)
   - Update static IP attachment

3. **Or Add Storage**
   - Attach block storage disk (additional cost)
   - Mount at `/var/www/nanosim/data`

## Security Recommendations

1. ✅ **Enable firewall** - Done in Step 3
2. ✅ **Use HTTPS** - Done in Step 8
3. ✅ **Keep system updated**: `sudo apt update && sudo apt upgrade`
4. ✅ **Change SSH port** (optional but recommended)
5. ✅ **Add basic auth** - See Optional section
6. ✅ **Regular backups** - Set up cron job for automatic backups
7. ✅ **Monitor logs** - Check for suspicious activity

## Next Steps

After deployment:

1. Test all functionality (viewing plots, notes, combining sims)
2. Upload initial simulation data
3. Set up automated backups
4. Monitor performance and costs
5. Consider adding monitoring (CloudWatch or free alternatives)
