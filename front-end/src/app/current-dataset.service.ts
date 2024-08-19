import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';

import { DatasetInfo } from './dataset-info';
import { DatasetInfoService } from './dataset-info.service';

@Injectable({
  providedIn: 'root'
})
export class CurrentDatasetService {

  private readonly currentDataset$: BehaviorSubject<DatasetInfo> = new BehaviorSubject<DatasetInfo> ({
    key: '',
    title: '',
    samples: [],
    integration_obs: []
  });
  readonly currentDataset: Observable<DatasetInfo> = this.currentDataset$.asObservable();

  private datasetInfo: DatasetInfo[] = [];

  constructor(private datasetInfoService: DatasetInfoService) {
    this.datasetInfoService.datasetInfo.subscribe(
      (res) => {
        this.datasetInfo = res;
        this.currentDataset$.next(this.datasetInfo[0]);
      },
    );
  }

  setCurrentDataset(newDataset: string): void {
    const currentDataset = this.datasetInfo.find(dataset => dataset.key == newDataset) || this.currentDataset$.getValue();
    this.currentDataset$.next(currentDataset);
  }
}
