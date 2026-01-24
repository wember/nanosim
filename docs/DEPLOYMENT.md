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
- [Step 7: Configure Nginx](#step-7-configure-nginx)
- [Step 8: Setup SSL with Let's Encrypt](#step-8-setup-ssl-with-lets-encrypt)
- [Step 9: Upload Simulation Data](#step-9-upload-simulation-data)
- [Step 10: Verify Deployment](#step-10-verify-deployment)
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
- Service: AWS Lightsail ($5/month)
- Server: Ubuntu 22.04 with Python/Flask
- Web Server: Nginx + Gunicorn
- SSL: Let's Encrypt (free, auto-renewal)
- Storage: 40GB SSD (included)

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
   - **Plan**: $5/month (1 GB RAM, 1 vCPU, 40 GB SSD)
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
workers = 2
worker_class = "sync"
timeout = 120
accesslog = "/var/log/nanosim/access.log"
errorlog = "/var/log/nanosim/error.log"
loglevel = "info"
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
```

### Start Gunicorn

```bash
sudo supervisorctl reread
sudo supervisorctl update
sudo supervisorctl start nanosim
sudo supervisorctl status nanosim  # Should show RUNNING
```

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

    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 300s;
        proxy_connect_timeout 75s;
    }

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
- Automatically update Nginx configuration
- Setup auto-renewal

### Verify SSL Auto-Renewal

```bash
# Test renewal
sudo certbot renew --dry-run

# Check renewal timer
sudo systemctl status certbot.timer
```

## Step 9: Upload Simulation Data

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

## Step 10: Verify Deployment

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

## Maintenance Tasks

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

## Cost Breakdown

- **Lightsail Instance**: $5.00/month (fixed)
- **Data Transfer**: Included (first 2 TB)
- **Static IP**: Included
- **SSL Certificate**: Free (Let's Encrypt)
- **DNS (Route 53)**: ~$0.50/month per hosted zone

**Total: ~$5-6/month**

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
