# pylint: disable=
"""
Log Utility Module
"""
# ============================== Standard Library Imports ==============================
import datetime
import time
from typing import List


# ============================== Custom Functions ==============================
def get_task_completion_message(start_time: float) -> List[str]:
    """
    Return the core task-completion message lines without borders.

    Parameters
    ----------
    start_time : float
        Timestamp at the beginning of the task (from time.time()).

    Returns
    -------
    List[str]
        List of plain log message lines (without borders).
    """
    end_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    duration = str(datetime.timedelta(seconds=round(time.time() - start_time)))

    return [
        "Task Finished Successfully.",
        f"End Time: {end_time} | Duration: {duration}",
    ]


def get_pipeline_completion_message(start_time: float) -> List[str]:
    """
    Return the core pipeline-completion message lines without borders.

    Parameters
    ----------
    start_time : float
        Timestamp at the beginning of the task (from time.time()).

    Returns
    -------
    List[str]
        List of plain log message lines (without borders).
    """
    end_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    duration = str(datetime.timedelta(seconds=round(time.time() - start_time)))

    return [
        "Pipeline Finished Successfully.",
        f"End Time: {end_time} | Duration: {duration}",
    ]
