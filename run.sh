#!/usr/bin/env bash
# ===========================================================================
# EpiNexus – single-command startup script
#
# Usage:
#   ./run.sh              Start all services (API + Streamlit frontend)
#   ./run.sh --api        Start API server only
#   ./run.sh --frontend   Start Streamlit frontend only
#   ./run.sh --stop       Stop all running EpiNexus services
# ===========================================================================

set -euo pipefail

# ── Project root ───────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Defaults (override via environment or .env) ───────────────────────────
API_HOST="${EPINEXUS_API_HOST:-127.0.0.1}"
API_PORT="${EPINEXUS_API_PORT:-8000}"
FRONTEND_PORT="${EPINEXUS_FRONTEND_PORT:-8501}"
LOG_DIR="${EPINEXUS_LOG_DIR:-$SCRIPT_DIR/logs}"
PID_DIR="${EPINEXUS_PID_DIR:-$SCRIPT_DIR/.pids}"
WORKERS="${EPINEXUS_WORKERS:-1}"

# ── Load .env if present ──────────────────────────────────────────────────
if [[ -f "$SCRIPT_DIR/.env" ]]; then
    set -a
    # shellcheck disable=SC1091
    source "$SCRIPT_DIR/.env"
    set +a
fi

# ── Helpers ───────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'  # No Color

info()  { echo -e "${GREEN}[INFO]${NC}  $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }

mkdir -p "$LOG_DIR" "$PID_DIR"

# ── Dependency checks ────────────────────────────────────────────────────
check_deps() {
    local missing=()
    command -v python3 &>/dev/null || missing+=("python3")
    python3 -c "import streamlit" 2>/dev/null || missing+=("streamlit (pip install streamlit)")
    python3 -c "import fastapi"  2>/dev/null || missing+=("fastapi  (pip install fastapi)")
    python3 -c "import uvicorn"  2>/dev/null || missing+=("uvicorn  (pip install uvicorn)")

    if [[ ${#missing[@]} -gt 0 ]]; then
        error "Missing dependencies:"
        for dep in "${missing[@]}"; do
            echo "  - $dep"
        done
        echo ""
        info "Install all dependencies with:"
        echo "  pip install -r requirements.txt"
        exit 1
    fi
}

# ── Service management ───────────────────────────────────────────────────
start_api() {
    if [[ -f "$PID_DIR/api.pid" ]] && kill -0 "$(cat "$PID_DIR/api.pid")" 2>/dev/null; then
        warn "API server already running (PID $(cat "$PID_DIR/api.pid"))"
        return
    fi

    info "Starting FastAPI server on http://${API_HOST}:${API_PORT} ..."
    python3 -m uvicorn app.main:app \
        --host "$API_HOST" \
        --port "$API_PORT" \
        --workers "$WORKERS" \
        --log-level info \
        >> "$LOG_DIR/api.log" 2>&1 &

    echo $! > "$PID_DIR/api.pid"
    info "API server started (PID $!)"
}

start_frontend() {
    if [[ -f "$PID_DIR/frontend.pid" ]] && kill -0 "$(cat "$PID_DIR/frontend.pid")" 2>/dev/null; then
        warn "Streamlit frontend already running (PID $(cat "$PID_DIR/frontend.pid"))"
        return
    fi

    info "Starting Streamlit frontend on http://localhost:${FRONTEND_PORT} ..."
    python3 -m streamlit run frontend/EpiNexus.py \
        --server.port "$FRONTEND_PORT" \
        --server.headless true \
        --browser.gatherUsageStats false \
        >> "$LOG_DIR/frontend.log" 2>&1 &

    echo $! > "$PID_DIR/frontend.pid"
    info "Streamlit frontend started (PID $!)"
}

stop_services() {
    local stopped=0
    for service in api frontend; do
        local pidfile="$PID_DIR/${service}.pid"
        if [[ -f "$pidfile" ]]; then
            local pid
            pid="$(cat "$pidfile")"
            if kill -0 "$pid" 2>/dev/null; then
                kill "$pid"
                info "Stopped $service (PID $pid)"
                stopped=$((stopped + 1))
            fi
            rm -f "$pidfile"
        fi
    done

    if [[ $stopped -eq 0 ]]; then
        info "No running EpiNexus services found."
    fi
}

wait_for_services() {
    info "─────────────────────────────────────────────"
    info "EpiNexus is starting up!"
    info ""
    info "  API:      http://${API_HOST}:${API_PORT}/docs"
    info "  Frontend: http://localhost:${FRONTEND_PORT}"
    info "  Logs:     $LOG_DIR/"
    info ""
    info "Press Ctrl+C to stop all services."
    info "─────────────────────────────────────────────"

    # Trap Ctrl-C to clean up
    trap 'echo ""; info "Shutting down..."; stop_services; exit 0' INT TERM

    # Wait for background processes
    wait
}

# ── Main ──────────────────────────────────────────────────────────────────
case "${1:-all}" in
    --api|-a)
        check_deps
        start_api
        wait_for_services
        ;;
    --frontend|-f)
        check_deps
        start_frontend
        wait_for_services
        ;;
    --stop|-s)
        stop_services
        ;;
    --help|-h)
        echo "Usage: $0 [--api | --frontend | --stop | --help]"
        echo ""
        echo "  (no args)    Start both API and Streamlit frontend"
        echo "  --api        Start API server only"
        echo "  --frontend   Start Streamlit frontend only"
        echo "  --stop       Stop all running services"
        ;;
    all|"")
        check_deps
        start_api
        sleep 1  # let API bind before frontend tries to connect
        start_frontend
        wait_for_services
        ;;
    *)
        error "Unknown option: $1"
        echo "Run '$0 --help' for usage."
        exit 1
        ;;
esac
