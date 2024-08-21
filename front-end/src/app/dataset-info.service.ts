import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';

import { BackendHttpClient } from './backend-http.client';
import { DatasetInfo } from './dataset-info';
import { DatasetInfoServiceInterface } from './dataset-info.service.interface';

@Injectable({
  providedIn: 'root'
})
export class DatasetInfoService implements DatasetInfoServiceInterface {

  private readonly datasetInfo$: BehaviorSubject<DatasetInfo[]> = new BehaviorSubject<DatasetInfo[]> ([]);
  readonly datasetInfo: Observable<DatasetInfo[]> = this.datasetInfo$.asObservable();

  constructor(private backendHttpClient: BackendHttpClient) { 
    this.setDatasetInfo();
  }

  /**
   * Asks the backendHttpClient to get the dataset info, and updates the dataset info when it receives a response.
   */
  async setDatasetInfo(): Promise<void> {
    const datasetInfo = await this.backendHttpClient.getDatasetInfo();
    this.datasetInfo$.next(datasetInfo);
  }
}
