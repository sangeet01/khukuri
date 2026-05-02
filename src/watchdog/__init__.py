"""
Khukuri Watchdog System — Continuous health and resistance monitoring.

Watchdogs:
1. SystemHealthWatchdog — Monitors system stability and module integrity.
2. ScientificHealthWatchdog — Monitors PINCER evolution and discovery quality.
3. DrugResistanceSentinel — Monitors deployed drugs against the resistance landscape.
"""

from .system_health import SystemHealthWatchdog, HealthReport, ComponentHealth
from .science_health import ScientificHealthWatchdog, ScienceHealthReport, MetricHealth
from .drug_sentinel import DrugResistanceSentinel, MonitoredDrug, ResistanceAlert

__all__ = [
    "SystemHealthWatchdog",
    "HealthReport",
    "ComponentHealth",
    "ScientificHealthWatchdog",
    "ScienceHealthReport",
    "MetricHealth",
    "DrugResistanceSentinel",
    "MonitoredDrug",
    "ResistanceAlert",
]
